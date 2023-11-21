
from all_functions import *

r = 0.065 # radius of neighborhood on unit sphere
coords = read_data_netcdf() #Read data from a netcdf file

coord_size = len(coords)

#Main Loop which iterates over all points and calls different functions
for id in range(coord_size):
    my_list = create_neighborhood(r, coords, id) #Myimplementation

    #print('mylist', my_list)
    #print(len(my_list))
    #print(coords[my_list[id]])

    xyz_r = getneighcoords(my_list, coords)
    #print("sub coords:", xyz_r)

    A, invA = constructA(xyz_r)

    #print('Matrix A:', A)

    cnd = linalg.cond(A)

    print('Condition number of A ', cnd)

    px, py, pz = projection(xyz_r)

    Dx, Dy, Dz = differentiation(invA,xyz_r,coords[id])


    print('***End of loop ', id, '***')



