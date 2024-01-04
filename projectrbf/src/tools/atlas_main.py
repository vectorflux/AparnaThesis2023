# Master Thesis Project

# This is the main python script for the project
# Step by step , it will call and execute all the functions.
#   1. Reads Coordinates list from Netcdf4 file
#   2. Provides coordinates input to Atlas Point Cloud script
#   3. Gets Indices back from Atlas Script
#   4. Creates A Matrix and its inverse
#   5. Creates Projection Matrix P : Px, Py, Pz
#   6. Creates Differentiation Matrix D : B * A inv
#   7. Time Stepping Loop
#   .
#   .
#   8. Run Test Cases
#   9. Publish Results

#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

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


atlas_start = time.time()
atlas.initialize() #initializes atlas and MPI



#Constants and Variables
myradius = 0.065 #can generalize later
allA = []
allinvA = []
nrj_size_list = []
allnearest = []
allDx = []
allDy = []
allDz = []
allpx = []
allpy = []
allpz = []
f = []


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

myfield = functionspace.create_field(name="swe_variables", variables=4, dtype = np.float64)
uvwh = atlas.make_view(myfield)
 

mycoords = functionspace.create_field(name="lonlat", variables=2, dtype = np.float64)
lonlat_sd = atlas.make_view(mycoords)
#lonlat_sd = lonlat

   
#print("Lonlat SD:", lonlat_sd)


### FINDING NEIGHBORS + INITIALIZATION
n_p = functionspace.size
search = Search(functionspace) #initializes the Search class with functionspace
xyz_np = np.zeros([n_p,3])

cnd = np.zeros(len(internal))

#loops over all the points in the subdomain
#each function will work on one j at a time

##################################333###############
#Big initialization loop
n=0
initloop_start = time.time()
for id in range(n_p):  # n_p
    if ghost[id] == 0:

        #nearest = find_neighbors(radius,functionspace) #nearest gives set of local indices
        nearest = search.nearest_indices_within_radius(id, myradius)
        allnearest = np.append(allnearest, nearest)
        nrj_size_list = np.append(nrj_size_list, len(nearest))
        xyz_r = getneighcoords(nearest, xyz) #Gets the list of neighborhood coordinates for further use
        #if(id == 0):
        #   print("nearest : ", nearest)


        #call function to create A, A inverse and save them
        A, invA = constructA(xyz_r)
        cnd[n] = linalg.cond(A)
        n=n+1
        allA = np.append(allA,A)
        allinvA = np.append(allinvA,invA)
        
        #xyz_np[id]=xyz[id]

        #call function to create D
        Xj = xyz[id]
        Dx, Dy, Dz = differentiation(invA,xyz_r,Xj)
        allDx = np.append(allDx, Dx)
        allDy = np.append(allDy, Dy)
        allDz = np.append(allDz, Dz)

        px, py, pz = projection(xyz_r)
        allpx = np.append(allpx, px)
        allpy = np.append(allpy, py)
        allpz = np.append(allpz, pz)

initloop_end = time.time()

#print("Time elapsed in the initialization loop): ", initloop_end-initloop_start, "seconds" )
# call function to  initialize fields
allD = np.column_stack([allDx, allDy, allDz])
allP = np.column_stack([allpx, allpy, allpz])

print("Condition number of A", np.mean(cnd))

#print("max of nrj size list : ", np.max(nrj_size_list))
#print("allnearest : ", allnearest)
#print("allinvA|size : " ,len(allinvA), "\n", allinvA)
#print("all Px |size:", len(allpx), "\n", allpx )
#print("all Py |size:", len(allpy), "\n", allpy )
#print("all Pz |size:", len(allpz), "\n", allpz )
#print("all Px : \n", allpx )
#print("all Py : \n", allpy )
#print("all Pz : \n", allpz )


#########################################################
def get_rk4_values(uvwh, f, dt, nrj_size_list, allnearest,xyz, allD, ghost):

    #Here uvwh has the size of entire subdomain n_p including Halo Region

    length = len(uvwh)
    rk_u = np.ndarray([length,4])
    rk_v = np.ndarray([length,4])
    rk_w = np.ndarray([length,4])
    rk_h = np.ndarray([length,4])
    uvwh0 = uvwh
    num = len(uvwh)
    #for i in range(len(uvwh)):
        ##f ghost[i]:
         #   uvwh0 = np.delete(uvwh, (i), axis = 0)
            #print("Halo region is present")

    #halo = [i for i in range(len(uvwh)) if ghost[i]==1]
    #uvwh0 = np.delete(uvwh)
    #uvwh0 = np.delete(uvwh, (halo), axis = 0)

    #print("New size of uvwh:", len(uvwh0))

    
    #d1 = dt*rhs(k1)
    rk_u[:,0], rk_v[:,0], rk_w[:,0], rk_h[:,0]  = construct_rhsd(nrj_size_list, allnearest, uvwh, f, xyz, allD,ghost)
    rk_u[:, 0] = rk_u[:, 0] * dt
    rk_v[:, 0] = rk_v[:, 0] * dt
    rk_w[:, 0] = rk_w[:, 0] * dt
    rk_h[:, 0] = rk_h[:, 0] * dt

    #******************  k2  **************#

    #k2 = uvwh0 + 0.5*d1
    #k = 0
    for i in range(num):
        if not ghost[i]:
            uvwh[i, 0] = uvwh0[i, 0] + (rk_u[i,0]/2)
            uvwh[i, 1] = uvwh0[i, 1] + (rk_v[i,0]/2)
            uvwh[i, 2] = uvwh0[i, 2] + (rk_w[i,0]/2)
            uvwh[i, 3] = uvwh0[i, 3] + (rk_h[i,0]/2)
            #k=k+1

        #Update values (Halo exchange)
    myfield.halo_dirty = True
    myfield.halo_exchange()

        # d2 = dt*rhs(k2)
    rk_u[:, 1], rk_v[:, 1], rk_w[:, 1], rk_h[:, 1] = construct_rhsd(nrj_size_list, allnearest, uvwh, f, xyz, allD, ghost)
    rk_u[:, 1] = rk_u[:, 1] * dt
    rk_v[:, 1] = rk_v[:, 1] * dt
    rk_w[:, 1] = rk_w[:, 1] * dt
    rk_h[:, 1] = rk_h[:, 1] * dt

    #*************   k3  *************#

    #k3 = uvwh0 +0.5*d2
    #k = 0
    for i in range(num):
        if not ghost[i]:

            uvwh[i, 0] = uvwh0[i, 0] + (rk_u[i, 1]/2)
            uvwh[i, 1] = uvwh0[i, 1] + (rk_v[i, 1]/2)
            uvwh[i, 2] = uvwh0[i, 2] + (rk_w[i, 1]/2)
            uvwh[i, 3] = uvwh0[i, 3] + (rk_h[i, 1]/2)
            #k=k+1

    #Update values (Halo exchange)
    myfield.halo_dirty = True
    myfield.halo_exchange()

     # d3 = dt*rhs(k3)
    rk_u[:, 2], rk_v[:, 2], rk_w[:, 2], rk_h[:, 2] = construct_rhsd(nrj_size_list, allnearest, uvwh, f, xyz, allD,ghost)
    rk_u[:, 2] = rk_u[:, 2] * dt
    rk_v[:, 2] = rk_v[:, 2] * dt
    rk_w[:, 2] = rk_w[:, 2] * dt
    rk_h[:, 2] = rk_h[:, 2] * dt

    # *************   k4  *************#

    #k4 = uvwh0 + d3
    #k = 0
    for i in range(num):
        if not ghost[i]:
            uvwh[i, 0] = uvwh0[i, 0] + (rk_u[i, 2])
            uvwh[i, 1] = uvwh0[i, 1] + (rk_v[i, 2])
            uvwh[i, 2] = uvwh0[i, 2] + (rk_w[i, 2])
            uvwh[i, 3] = uvwh0[i, 3] + (rk_h[i, 2])
            #k=k+1

    #Update values (Halo exchange)
    myfield.halo_dirty = True 
    myfield.halo_exchange()

    # calculation of rhs with k4 as input
    rk_u[:, 3], rk_v[:, 3], rk_w[:, 3], rk_h[:, 3] = construct_rhsd(nrj_size_list, allnearest, uvwh, f, xyz, allD, ghost)

    # d4 = dt*rhs(k4)
    rk_u[:, 3] = rk_u[:, 3] * dt
    rk_v[:, 3] = rk_v[:, 3] * dt
    rk_w[:, 3] = rk_w[:, 3] * dt
    rk_h[:, 3] = rk_h[:, 3] * dt

    return rk_u, rk_v, rk_w, rk_h
###################################################################

#function to initialize fields
uvwh, f = set_initial_conditions(uvwh, xyz, n_p, ghost,lonlat)

for n in range(n_p):
  if not ghost[n]:
    lonlat_sd[n,0] = lonlat[n,0]
    lonlat_sd[n,1] = lonlat[n,1]
    
  else:
     lonlat_sd[n,0] = lonlat[n,1] = 0.0

#print("Range: ", np.min)



timeloop_start = time.time()

n_timesteps = 500 #12 days (12 * 86400/1200)

for i in range(n_timesteps):

    #Update values (Halo exchange)
    myfield.halo_dirty = True  # 
    myfield.halo_exchange()

    dt = 2 #2 mins
    #halo = [i for i in range(len(uvwh)) if ghost[i]==1]
    
    uvwh0 = uvwh #with halo region


    # function to get arrays of rk values for u, v, w, h
    #rk_u will be the size of np points

    rk_u, rk_v, rk_w, rk_h = get_rk4_values(uvwh, f, dt, nrj_size_list, allnearest,xyz, allD, ghost)
    
    #print("rk4 in the time loop",i ,"are: ", rk_u, rk_v, rk_w, rk_h)
    
    #halo = [i for i in range(len(uvwh)) if ghost[i]==1]
    #uvwh0 = uvwh = np.delete(uvwh, (halo), axis = 0)
    #l = 0 
    for k in range(n_p):
        if not ghost[k]:
            
            #Calculate uvwh at next timestep (overwriting): need all RHS values
            uvwh[k,0] = uvwh0[k,0] + (rk_u[k,0] + rk_u[k,3])/6 + (rk_u[k,1]+rk_u[k,2])/3
            uvwh[k,1] = uvwh0[k,1] + (rk_v[k,0] + rk_v[k,3])/6 + (rk_v[k,1]+rk_v[k,2])/3
            uvwh[k,2] = uvwh0[k,2] + (rk_w[k,0] + rk_w[k,3])/6 + (rk_w[k,1]+rk_w[k,2])/3
            uvwh[k,3] = uvwh0[k,3] + (rk_h[k,0] + rk_h[k,3])/6 + (rk_h[k,1]+rk_h[k,2])/3

            #l =l+1

    #print("uvwh in the time loop",i ,"are: ", uvwh)

    #print("Final value of k:", k,"for loop of i", i)

    #uvwh[:,1] = uvwh0[:,1] + (rk_v[:,0] + rk_v[:,3])/6 + (rk_v[:,1]+rk_v[:,2])/3
    #update uvwh with new values
    #uvwh = uvwh_n


timeloop_end = time.time()
#print("Time taken by the time loop: ", timeloop_end-timeloop_start, "seconds" )
#print("rk4 final are: ", rk_u, rk_v, rk_w, rk_h)

#print("Final values of uvwh: " , uvwh)
#print("Range of u: ", np.min(uvwh[:,0]), np.max(uvwh[:,0]))
#print("Range of v: ", np.min(uvwh[:,1]), np.max(uvwh[:,1]))
#print("Range of w: ", np.min(uvwh[:,2]), np.max(uvwh[:,2]))
#print("Range of h: ", np.min(uvwh[:,3]), np.max(uvwh[:,3]))

myfield.halo_dirty = True
myfield.halo_exchange()

mycoords.halo_dirty = True
mycoords.halo_exchange()

print("Final values of uvwh: " , uvwh)

fs = myfield.functionspace
field_global = fs.create_field_global(name="swe_variables_global", variables=4,dtype=np.float64)
functionspace.gather(myfield, field_global)
uvwh_global = atlas.make_view(field_global)

fs = mycoords.functionspace
lonlat_global = fs.create_field_global(name="lonlat_global", variables=2)
functionspace.gather(mycoords, lonlat_global)
lonlat_global_coords = atlas.make_view(lonlat_global)

if functionspace.part ==0:
    #print("Lonlat Global:", lonlat_global_coords)
    #print("UVWH Global:", uvwh_global)
    plot_global(uvwh_global,lonlat_global_coords)



atlas.finalize()

atlas_end = time.time()

print("Total Time elapsed (Atlas-finalize): ", atlas_end-atlas_start)
