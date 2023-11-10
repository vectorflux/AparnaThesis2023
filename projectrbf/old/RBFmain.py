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

#Authors: Aparna Devulapalli, William Sawyer, Willem Deconinck


from projectrbf.src.tools.operator_matrices import *
from projectrbf.src.tools.read_netcdf_file import *
from projectrbf.src.atlas.atlas_func_interface import *


def execute_main():



    #Define radius as per the size of neighborhood required
    r = 0.065

    #Get data from netcdf in both dgrees and cartesian forms
    lonlat, coords = read_data_netcdf()

    #Get nearest from Atlas for all indices
    #need to modify
    indices_list = find_neighbors(r, lonlat)

    #Get inverse of A
    invA = constructa(indices_list,coords)

    #create the differentiator matrices
    differentiation(sub_coords, invA, coords[id])




    #print('time elapsed:', end-start)
#print(A)

if __name__ == '__main__':
    print('Running RBF Main')
    execute_main()
