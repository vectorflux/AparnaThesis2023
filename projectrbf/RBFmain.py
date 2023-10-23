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


import time
import numpy as np
from rbf_matrices import *
from operator_matrices import *
from read_netcdf_file import *
import numpy.linalg as linalg


# Use Atlas in place of Create Neighborhoods
def create_neighborhood(r, coords, my_i):
    coord_size = len(coords)
    #print('coordsize ', coord_size)
    #print('r', r)
    list = []
    #ctrlist = []

    #for i in range(coord_size):
    ctr = 0

    for j in range(coord_size):
        if (j != my_i):
            radius = eucl_norm(coords[my_i],coords[j])
            #print('entering first if')
            if (radius<=r):
                #print('entering second if')
                ctr+=1
                list.append(j)

    return list


def execute_main():
    start = time.time()
    coords = read_data_netcdf()
    #print(coords[:][:])
    #r = 400000
    r = 0.065
    coord_size = len(coords)


    #Main Loop which iterates over all points and calls different functions
    for id in range(coord_size):

        my_list = create_neighborhood(r, coords, id) #Myimplementation
        #my_list = find_neighbors(r) #atlasimplemntation

        #print('mylist', my_list)
        #print(len(my_list))
        #print(coords[my_list[id]])

        sub_coords = getsubdomain(my_list, coords)
        #print("sub coords:", sub_coords)

        A = constructa(sub_coords)
        #print('Matrix A:', A)

        #cnd = linalg.cond(A)
        #print('Condition number of A ', cnd)

        InvA = inverta(A)
        #p = projection(coords,my_list)

        differentiation(sub_coords,InvA,coords[id])


        #print('P Matrix:', p)
        print('***End of loop ', id, '***')

        end = time.time()

    print('time elapsed:', end-start)
#print(A)

if __name__ == '__main__':
    print('Running read points buildA file')
    execute_main()
