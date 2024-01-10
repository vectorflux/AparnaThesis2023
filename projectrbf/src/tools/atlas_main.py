# Master Thesis Project

# This is the main python script for the project
# Step by step , it will call and execute all the functions.
#   1. Reads Coordinates list from Netcdf4 file
#   2. Provides coordinates input to Atlas Point Cloud script
#   3. Gets Indices back from Atlas Script
#   4. Creates A Matrix and its inverse
#   5. Creates Differentiation Matrix D : B * A inv
#   6. Initializes fields with test case conditions
#   7. RK4 Time Stepping Loop
#   8. Run Test Cases
#   9. Publish Results


#!/usr/bin/env python3
import atlas4py as atlas
import numpy as np
import time

from operator_matrices import *
from atlas_func_interface import *
from read_netcdf_file import *
from A_matrices import *
from initializefields import *
from visualization import *
from rk4 import *
from construct_rhsd import *


atlas_start = time.time()
atlas.initialize() #initializes atlas and MPI

#Constants and Variables
myradius = 0.065 #can generalize later
nrj_size_list = []
allnearest = []
allDx = []
allDy = []
allDz = []


lonlat = read_data_netcdf() #Read data from a netcdf file
grid = atlas.UnstructuredGrid(lonlat[:, 0], lonlat[:, 1]) #Create Unstructured Grid

#Create functionspace + partitioning #levels in the pointcloud implementation
functionspace = atlas.functionspace.PointCloud(grid, halo_radius=myradius * 2, geometry="UnitSphere")

###Access fields as numpy arrays, if needed
lonlat = atlas.make_view(functionspace.lonlat) # longitude latitude in degrees
ghost = atlas.make_view(functionspace.ghost) # ghost: 1 for ghost, 0 for owned
partition = atlas.make_view(functionspace.partition) #partition owning point (0-based)
remote_index = atlas.make_view(functionspace.remote_index) # local index on partition owning point (careful, 1-based when atlas_HAVE_FORTRAN)
global_index = atlas.make_view(functionspace.global_index) # global index across partitions (always 1-based)



xyz = getcartesian(lonlat)
###Create Fields
#first dimension - size of the functionspace (rows)
#second dimension - number of variables (columns)
internal = [ i for i in range(functionspace.size) if ghost[i] == 0 ]
#field = functionspace.create_field(name="mycoords",variables=3,levels = 10, dtype = np.float64)

myfield = functionspace.create_field(name="swe_variables", variables=5, dtype = np.float64)
uvwh = atlas.make_view(myfield)

mycoords = functionspace.create_field(name="lonlat", variables=2, dtype = np.float64)
lonlat_sd = atlas.make_view(mycoords)

n_p = functionspace.size
search = Search(functionspace) #initializes the Search class with functionspace

initloop_start = time.time()
###################  Loop to construct RBF Operator matrice  ############################
initloop_start = time.time()
for index in range(n_p):  # n_p
    if ghost[index] == 0:
        
        nearest = search.nearest_indices_within_radius(index, myradius)
        allnearest = np.append(allnearest, nearest)
        nrj_size_list = np.append(nrj_size_list, len(nearest))
        xyz_r = getneighcoords(nearest, xyz) #Gets the list of neighborhood cartesian coordinates for further use

        #call function to create A, A inverse and save them
        A, invA = constructA(xyz_r)
        
        Xj = xyz[index]
        #call function to create gradient operators
        Dx, Dy, Dz = differentiation(invA,xyz_r,Xj)
        
        allDx = np.append(allDx, Dx)
        allDy = np.append(allDy, Dy)
        allDz = np.append(allDz, Dz)

initloop_end = time.time()
allD = np.column_stack([allDx, allDy, allDz])


###################################################################

#function call to get initial conditions
fields = set_initial_conditions(xyz, n_p, ghost,lonlat)

#fill values in the functionspace variables
for n in range(n_p):
  if not ghost[n]:
    lonlat_sd[n,0] = lonlat[n,0]
    lonlat_sd[n,1] = lonlat[n,1]
    
    uvwh[n,0] = fields[n,0]
    uvwh[n,1] = fields[n,1]
    uvwh[n,2] = fields[n,2]
    uvwh[n,3] = fields[n,3]
    uvwh[n,4] = fields[n,4]
     
  else:
     lonlat_sd[n,0] = lonlat[n,1] = 0.0
     uvwh[n,0] = 0.0
     uvwh[n,1] = 0.0
     uvwh[n,2] = 0.0
     uvwh[n,3] = 0.0
     uvwh[n,4] = 0.0


#Update values (Halo exchange)
myfield.halo_dirty = True
myfield.halo_exchange()

mycoords.halo_dirty = True
mycoords.halo_exchange()

timeloop_start = time.time()

n_timesteps = 1 #12 days (12 * 86400/1200)


#########################################################
def get_rk4_values(uvwh, dt, nrj_size_list, allnearest,xyz, allD, ghost):

    #Here uvwh has the size of entire subdomain n_p including Halo Region
    num = len(uvwh)

    rk_u = np.zeros((num,4))
    rk_v = np.zeros((num,4))
    rk_w = np.zeros((num,4))
    rk_h = np.zeros((num,4))
    uvwh0 = uvwh
    num = len(uvwh)

    #d1 = dt*rhs(k1)
    rk_u[:,0], rk_v[:,0], rk_w[:,0], rk_h[:,0]  = construct_rhsd(nrj_size_list, allnearest, uvwh, xyz, allD, ghost)
    rk_u[:, 0] = rk_u[:, 0] * dt
    rk_v[:, 0] = rk_v[:, 0] * dt
    rk_w[:, 0] = rk_w[:, 0] * dt
    rk_h[:, 0] = rk_h[:, 0] * dt


    #******************  k2  **************#

    #k2 = uvwh0 + 0.5*d1
    for i in range(num):
        uvwh[i, 0] = uvwh0[i, 0] + (rk_u[i,0]/2)
        uvwh[i, 1] = uvwh0[i, 1] + (rk_v[i,0]/2)
        uvwh[i, 2] = uvwh0[i, 2] + (rk_w[i,0]/2)
        uvwh[i, 3] = uvwh0[i, 3] + (rk_h[i,0]/2)


    #Update values (Halo exchange)
    myfield.halo_dirty = True
    myfield.halo_exchange()

    # d2 = dt*rhs(k2)
    rk_u[:, 1], rk_v[:, 1], rk_w[:, 1], rk_h[:, 1] = construct_rhsd(nrj_size_list, allnearest, uvwh, xyz, allD, ghost)
    rk_u[:, 1] = rk_u[:, 1] * dt
    rk_v[:, 1] = rk_v[:, 1] * dt
    rk_w[:, 1] = rk_w[:, 1] * dt
    rk_h[:, 1] = rk_h[:, 1] * dt

    #*************   k3  *************#

    #k3 = uvwh0 +0.5*d2
    for i in range(num):
        uvwh[i, 0] = uvwh0[i, 0] + (rk_u[i, 1]/2)
        uvwh[i, 1] = uvwh0[i, 1] + (rk_v[i, 1]/2)
        uvwh[i, 2] = uvwh0[i, 2] + (rk_w[i, 1]/2)
        uvwh[i, 3] = uvwh0[i, 3] + (rk_h[i, 1]/2)

    #Update values (Halo exchange)
    myfield.halo_dirty = True
    myfield.halo_exchange()

     # d3 = dt*rhs(k3)
    rk_u[:, 2], rk_v[:, 2], rk_w[:, 2], rk_h[:, 2] = construct_rhsd(nrj_size_list, allnearest, uvwh, xyz, allD, ghost)
    rk_u[:, 2] = rk_u[:, 2] * dt
    rk_v[:, 2] = rk_v[:, 2] * dt
    rk_w[:, 2] = rk_w[:, 2] * dt
    rk_h[:, 2] = rk_h[:, 2] * dt

    # *************   k4  *************#

    #k4 = uvwh0 + d3
    for i in range(num):
        uvwh[i, 0] = uvwh0[i, 0] + (rk_u[i, 2])
        uvwh[i, 1] = uvwh0[i, 1] + (rk_v[i, 2])
        uvwh[i, 2] = uvwh0[i, 2] + (rk_w[i, 2])
        uvwh[i, 3] = uvwh0[i, 3] + (rk_h[i, 2])

    #Update values (Halo exchange)
    myfield.halo_dirty = True
    myfield.halo_exchange()

    # calculation of rhs with k4 as input
    rk_u[:, 3], rk_v[:, 3], rk_w[:, 3], rk_h[:, 3] = construct_rhsd(nrj_size_list, allnearest, uvwh, xyz, allD, ghost)

    # d4 = dt*rhs(k4)
    rk_u[:, 3] = rk_u[:, 3] * dt
    rk_v[:, 3] = rk_v[:, 3] * dt
    rk_w[:, 3] = rk_w[:, 3] * dt
    rk_h[:, 3] = rk_h[:, 3] * dt

    return rk_u, rk_v, rk_w, rk_h


####################################################################
##### TIME STEPPING LOOP #######
####################################################################
for i in range(n_timesteps):

    rk_u = np.zeros((n_p,4))
    rk_v = np.zeros((n_p,4))
    rk_w = np.zeros((n_p,4))
    rk_h = np.zeros((n_p,4))
    
    dt = 1200 #20 mins
    
    uvwh0 = uvwh #store previous time step values

    rk_u, rk_v, rk_w, rk_h = get_rk4_values(uvwh, dt, nrj_size_list, allnearest,xyz, allD, ghost)
    
    for k in range(n_p):
        if not ghost[k]:
        
            #Calculate uvwh at next timestep (overwriting): need all RHS values
            uvwh[k,0] = uvwh0[k,0] + ((rk_u[k,0] + rk_u[k,3]) + 2*(rk_u[k,1]+rk_u[k,2]))/6
            uvwh[k,1] = uvwh0[k,1] + ((rk_v[k,0] + rk_v[k,3]) + 2*(rk_v[k,1]+rk_v[k,2]))/6
            uvwh[k,2] = uvwh0[k,2] + ((rk_w[k,0] + rk_w[k,3]) + 2*(rk_w[k,1]+rk_w[k,2]))/6
            uvwh[k,3] = uvwh0[k,3] + ((rk_h[k,0] + rk_h[k,3]) + 2*(rk_h[k,1]+rk_h[k,2]))/6

    #Update values after each time step (Halo exchange)
    myfield.halo_dirty = True
    myfield.halo_exchange()

timeloop_end = time.time()


#Update values (Halo exchange)
myfield.halo_dirty = True
myfield.halo_exchange()

mycoords.halo_dirty = True
mycoords.halo_exchange()


#create global variables and gather individual process data to the global variable at process 0
fs = myfield.functionspace
field_global = fs.create_field_global(name="swe_variables_global", variables=5,dtype=np.float64)
functionspace.gather(myfield, field_global)
uvwh_global = atlas.make_view(field_global)

mc = mycoords.functionspace
lonlat_global = mc.create_field_global(variables=2)
functionspace.gather(mycoords, lonlat_global)
lonlat_global_coords = atlas.make_view(lonlat_global)


if functionspace.part ==0:
    plot_global(uvwh_global,lonlat_global_coords)

print("Fin.")
atlas.finalize()

