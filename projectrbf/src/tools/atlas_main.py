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
#from rk4 import *
from construct_rhsd_vector import *


atlas_start = time.time()
atlas.initialize() #initializes atlas and MPI

#Constants and Variables
myradius = 0.075#can generalize later
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

myfield = functionspace.create_field(name="swe_variables", variables=4, dtype = np.float64)  # fifth variable is Coriolis force
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
        A, invA = constructA(xyz_r,myradius)
        
        Xj = xyz[index]
        #call function to create gradient operators
        Dx, Dy, Dz = differentiation(invA,xyz_r,Xj,myradius)
        
        allDx = np.append(allDx, Dx)
        allDy = np.append(allDy, Dy)
        allDz = np.append(allDz, Dz)

initloop_end = time.time()
allD = np.column_stack([allDx, allDy, allDz])


###################################################################

#function call to get initial conditions
fields, f = set_initial_conditions(xyz, n_p, ghost,lonlat)

#fill values in the functionspace variables
for n in range(n_p):
  if not ghost[n]:
    lonlat_sd[n,0] = lonlat[n,0]
    lonlat_sd[n,1] = lonlat[n,1]
    
    uvwh[n,0] = fields[n,0]
    uvwh[n,1] = fields[n,1]
    uvwh[n,2] = fields[n,2]
    uvwh[n,3] = fields[n,3]
     
  else:
     lonlat_sd[n,0] = lonlat[n,1] = 0.0
     uvwh[n,0] = 0.0
     uvwh[n,1] = 0.0
     uvwh[n,2] = 0.0
     uvwh[n,3] = 0.0

plot_global( uvwh, lonlat )   # plot initial conditions

#Update values (Halo exchange)
myfield.halo_dirty = True
myfield.halo_exchange()

mycoords.halo_dirty = True
mycoords.halo_exchange()

timeloop_start = time.time()


#########################################################
def get_rk4_values(uvwh, dt, nrj_size_list, allnearest,xyz, allD, ghost):

    #Here uvwh has the size of entire subdomain n_p including Halo Region
    num = len(uvwh)

    d1   = uvwh
    d2   = uvwh
    d3   = uvwh
    d4   = uvwh

    #d1 = dt*rhs(k1)
    d1  = dt * construct_rhsd_vector(nrj_size_list, allnearest, myradius, uvwh, f, xyz, allD, ghost) 

    #******************  k2  **************#

    #Update values (Halo exchange)
    myfield.halo_dirty = True
    myfield.halo_exchange()

    # k2 = uvwh0 + 0.5*d1
    # d2 = dt*rhs(k2)
    d2 = dt * construct_rhsd_vector(nrj_size_list, allnearest, myradius, (uvwh+d1/2), f, xyz, allD, ghost) 

    #*************   k3  *************#


    #Update values (Halo exchange)
    myfield.halo_dirty = True
    myfield.halo_exchange()

    # k3 = uvwh0 +0.5*d2
    # d3 = dt*rhs(k3)
    d3 = dt * construct_rhsd_vector(nrj_size_list, allnearest, myradius, (uvwh+d2/2), f, xyz, allD, ghost)

    # *************   k4  *************#

    #Update values (Halo exchange)
    myfield.halo_dirty = True
    myfield.halo_exchange()

    # calculation of rhs with k4 as input

    # k4 = uvwh0 + d3
    # d4 = dt*rhs(k4)
    d4 = dt * construct_rhsd_vector(nrj_size_list, allnearest, myradius, (uvwh+d3), f, xyz, allD, ghost)

    return d1, d2, d3, d4


####################################################################
##### TIME STEPPING LOOP #######
####################################################################

n_timesteps = 1 #12 days (12 * 86400/dt)

for i in range(n_timesteps):

    dt = 900 # 15 mins
    
    d1, d2, d3, d4 = get_rk4_values(uvwh, dt, nrj_size_list, allnearest,xyz, allD, ghost)
    
    drhs = ((d1 + d4) + 2*(d2+d3))/6
    plot_global(drhs,lonlat)
    for k in range(n_p):
        if not ghost[k]:
        
            #Calculate uvwh at next timestep (overwriting): need all RHS values
            uvwh[k,:] += drhs[k,:]

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
    plot_global( uvwh, lonlat )   # plot initial conditions

#    plot_global(uvwh_global,lonlat_global_coords)
    #gidx = atlas.make_view(functionspace.global_index)
#lonlat = atlas.make_view(functionspace.lonlat)if functionspace.part ==0:    for i in range(functionspace.size):
    #if not ghost[i]:
        #hlocal = uvwh[i,3]
        #hglobal = uvwh_global[gidx[i]-1,3]
        #if not (hlocal == hglobal):
            #print(hlocal, hglobal) 


print("Fin.")
atlas.finalize()

