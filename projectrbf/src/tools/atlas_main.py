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

start = time.time()
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
#first dimension - size of the functionspace
#second dimension - number of variables

#field = functionspace.create_field(name="mycoords",variables=3,levels = 10, dtype = np.float64)

myfield = functionspace.create_field(name="swe_variables", variables=4, dtype = np.float64)
uvwh = atlas.make_view(myfield)
#u = uvwh[:,0]
#v = uvwh[:,1]
#w = uvwh[:,2]
#h = uvwh[:,3]

### FINDING NEIGHBORS + INITIALIZATION

n_p = functionspace.size
search = Search(functionspace) #initializes the Search class with functionspace


#Ru = np.zeros(len(n_p))
#Rv = np.zeros(len(n_p))
#Rw = np.zeros(len(n_p))
#Rh = np.zeros(len(n_p))


#loops over all the points in the subdomain
#each function will work on one j at a time
#Big initialization loop

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



# call function to  initialize fields
allD = np.column_stack([allDx, allDy, allDz])
allP = np.column_stack([allpx, allpy, allpz])

#print("max of nrj size list : ", np.max(nrj_size_list))
#print("allDs are : ", allD)
#print("allinvA|size : " ,len(allinvA), "\n", allinvA)
#print("all Px |size:", len(allpx), "\n", allpx )
#print("all Py |size:", len(allpy), "\n", allpy )
#print("all Pz |size:", len(allpz), "\n", allpz )

#print("all Px : \n", allpx )
#print("all Py : \n", allpy )
#print("all Pz : \n", allpz )
#Ru, Rv, Rw, Rh = construct_rhs(uvwh,allD,xyz, nrj_size_list)



#Get all As
#Get all Ds




    # stage 1 : Get all A inverses and all D matrices




#time loop

    #iterate over subdomain
# stage 2 :




#########################################################



#function to initialize fields
uvwh = set_initial_conditions(uvwh, xyz, n_p, ghost)

myfield.halo_dirty = True  # if set to False, following halo_exchange will have no effect. Note it was already True upon create_field
myfield.halo_exchange()
vt = validate_halo_exchange(uvwh, xyz, n_p,ghost)
print("Halo exchange passed with vt = ", vt)

print("uvwh values are:\n" , uvwh )


#function to create rhsd matrix


#function to create rhs for u,v,w,h

#Ru, Rv, Rw, Rh = construct_rhs(uvwh,allD, allP,xyz, nrj_size_list)


#Time loop

#for i in range(n_timesteps):

    #Update values (Halo exchange)




    #Calculate uvwh at next timestep : need all RHS values




atlas.finalize()

end = time.time()

print("Time elapsed: ", end-start, "seconds" )



