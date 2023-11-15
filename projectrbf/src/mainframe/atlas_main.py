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


import atlas4py as atlas
import numpy as np
from MyThesis2023.projectrbf.src.tools.A_matrices import *
from MyThesis2023.projectrbf.src.tools.operator_matrices import *
from MyThesis2023.projectrbf.src.tools.atlas_func_interface import *
from MyThesis2023.projectrbf.src.tools.read_netcdf_file import *


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
u = uvwh[:,0]
v = uvwh[:,1]
w = uvwh[:,2]
h = uvwh[:,3]

### FINDING NEIGHBORS + INITIALIZATION

n_p = functionspace.size
search = Search(functionspace) #initializes the Search class with functionspace

#loops over all the points in the subdomain
#each function will work on one j at a time
#Big initialization loop

for id in range(n_p):  # n_p
    if ghost[id] == 0:

        #nearest = find_neighbors(radius,functionspace) #nearest gives set of local indices
        nearest = search.nearest_indices_within_radius(id, myradius)
        allnearest[id] = nearest
        nrj_size_list = np.append(nrj_size_list, len(nearest))
        xyz_r = getneighcoords(nearest, xyz) #Gets the list of neighborhood coordinates for further use


        #call function to create A, A inverse and save them
        A, invA = constructA(xyz_r)
        allA = np.append(allA,A)
        allinvA = np.append(allinvA,invA)


        #call function to create D
        Xj = xyz[id]
        Dx, Dy, Dz = differentiation(invA,xyz_r,Xj)
        allDx = np.append(allDx,Dx)
        allDy = np.append(allDy, Dy)
        allDz = np.append(allDz, Dz)

        px, py, pz = projection(xyz_r)
        allpx = np.append(allpx, px)
        allpy = np.append(allpy, py)
        allpz = np.append(allpz, pz)

        # call function to  initialize fields

allD = np.column_stack(allDx, allDy, allDz)

Ru, Rv, Rw, Rh = construct_rhs(uvwh,allD,xyz, nrj_size_list)



#Get all As
#Get all Ds




    # stage 1 : Get all A inverses and all D matrices




#time loop

    #iterate over subdomain
# stage 2 :




#########################################################






atlas.finalize()




