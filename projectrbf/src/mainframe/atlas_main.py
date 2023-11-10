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
import math
import numpy as np
from projectrbf.src.tools import *
from projectrbf.src.atlas.atlas_func_interface import *


atlas.initialize() #initializes atlas and MPI

resolution = 0.065

#Read data from a netcdf file
xyz, lonlat = read_data_netcdf()

#Create Unstructured Grid
grid = atlas.UnstructuredGrid(lonlat[:, 0], lonlat[:, 1])

#Create functionspace + partitioning
functionspace = atlas.functionspace.PointCloud(grid, halo_radius=resolution * 2, geometry="UnitSphere")


###Access fields as numpy arrays, if needed
lonlat = atlas.make_view(functionspace.lonlat) # longitude latitude in degrees
ghost = atlas.make_view(functionspace.ghost) # ghost: 1 for ghost, 0 for owned
partition = atlas.make_view(functionspace.partition) #partition owning point (0-based)
remote_index = atlas.make_view(functionspace.remote_index) # local index on partition owning point (careful, 1-based when atlas_HAVE_FORTRAN)
global_index = atlas.make_view(functionspace.global_index) # global index across partitions (always 1-based)

field = functionspace.create_field(name="mycoords")
xyz = atlas.make_view(field)
xyz = getcartesian(lonlat)

### FINDING NEIGHBORS + INITIALIZATION

allnearest = np.zeros([n_p, 50])
n_p = functionspace.size
search = Search(functionspace) #initializes the Search class with functionspace
radius = resolution




#loops over all the points in the subdomain
for id in range(n_p):  # n_p
    if ghost[id] == 0:


        nearest = find_neighbors(radius,functionspace) #nearest gives set of local indices
        allnearest[id] = nearest



        initpoints(nearest,xyz)



#########################################################






atlas.finalize()




