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

#field = functionspace.create_field(name="mycoords",variables=3,levels = 10, dtype = np.float64)

myfield = functionspace.create_field(name="swe_variables", variables=4, dtype = np.float64)
uvwh = atlas.make_view(myfield)

### FINDING NEIGHBORS + INITIALIZATION

n_p = functionspace.size
search = Search(functionspace) #initializes the Search class with functionspace


#loops over all the points in the subdomain
#each function will work on one j at a time

##################################333###############
#Big initialization loop

initloop_start = time.time()
for id in range(n_p):  # n_p
    if ghost[id] == 0:

        #nearest = find_neighbors(radius,functionspace) #nearest gives set of local indices
        nearest = search.nearest_indices_within_radius(id, myradius)
        allnearest = np.append(allnearest, nearest)
        nrj_size_list = np.append(nrj_size_list, len(nearest))
        xyz_r = getneighcoords(nearest, xyz) #Gets the list of neighborhood coordinates for further use


        #call function to create A, A inverse and save them
        A, invA = constructA(xyz_r)
        allA = np.append(allA,A)
        allinvA = np.append(allinvA,invA)


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

print("Time elapsed in the initialization loop): ", initloop_end-initloop_start, "seconds" )
# call function to  initialize fields
allD = np.column_stack([allDx, allDy, allDz])
allP = np.column_stack([allpx, allpy, allpz])



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

#myfield.halo_dirty = True  # if set to False, following halo_exchange will have no effect. Note it was already True upon create_field
#myfield.halo_exchange()
#vt = validate_halo_exchange(uvwh, xyz, n_p,ghost)
#print("Halo exchange passed with vt = ", vt)
#print("uvwh values are:\n" , uvwh )

#function to create rhs for u,v,w,h
#Ru, Rv, Rw, Rh = construct_rhsd(nrj_size_list, allnearest, uvwh, xyz, allD)
#print("Ru:",Ru, "\nRv:",Rv, "\nRw:",Rw, "\nRh:",Rh )
#print("******size of Ru: ", len(Ru))

#function to initialize fields
uvwh = set_initial_conditions(uvwh, xyz, n_p, ghost)

#Time loop

timeloop_start = time.time()

t = 0
tot_t = 100
dt = 0.1
n_timesteps = tot_t/dt

for i in range(n_timesteps):

    #Update values (Halo exchange)
    myfield.halo_dirty = True  # if set to False, following halo_exchange will have no effect. Note it was already True upon create_field
    myfield.halo_exchange()

    dt = 0.1
    uvwh0 = uvwh


    # function to get arrays of rk values for u, v, w, h
    rk_u, rk_v, rk_w, rk_h = get_rk4_values(uvwh, dt, nrj_size_list, allnearest,xyz, allD)

    #Calculate uvwh at next timestep (overwriting): need all RHS values
    uvwh[:,0] = uvwh0[:,0] + (rk_u[:,0] + rk_u[:,3])/6 + (rk_u[:,1]+rk_u[:,2])/3
    uvwh[:,1] = uvwh0[:,1] + (rk_v[:,0] + rk_v[:,3])/6 + (rk_v[:,1]+rk_v[:,2])/3
    uvwh[:,2] = uvwh0[:,2] + (rk_w[:,0] + rk_w[:,3])/6 + (rk_w[:,1]+rk_w[:,2])/3
    uvwh[:,3] = uvwh0[:,3] + (rk_h[:,0] + rk_h[:,3])/6 + (rk_h[:,1]+rk_h[:,2])/3

    #update uvwh with new values
    #uvwh = uvwh_n


timeloop_end = time.time()
print("Time taken by the time loop: ", timeloop_end-timeloop_start, "seconds" )

print("Final values of uvwh: " , uvwh)

atlas.finalize()

atlas_end = time.time()

print("Total Time elapsed (Atlas-finalize): ", atlas_end-atlas_start, "seconds" )



