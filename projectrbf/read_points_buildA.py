#Consider writing indices to a text file or have as many number of variables



import netCDF4 as nc
import time
import numpy as np
from rbf_matrices import *
from operator_matrices import *
import numpy.linalg as linalg

def read_data_netcdf():

    f = nc.Dataset('/home/dlaparna/thesis/projectrbf/GRID_FILES/grid.nc')
    #print(f.variables.keys())


    clon = f.variables['clon']
    clat = f.variables['clat']
    lon = clon[:]
    lat = clat[:]
    #print('number of longitudes', len(lon))
    #print('number of latitudes', len(lat))

    cx = f.variables['cartesian_x_vertices']
    cy = f.variables['cartesian_y_vertices']
    cz = f.variables['cartesian_z_vertices']


    R = 1 # radius of the earth in meters
    x = R * np.cos(lat) * np.cos(lon)
    y = R * np.cos(lat) * np.sin(lon)
    z = R * np.sin(lat)

    #print(len(x))

    coordxy = np.column_stack([x,y,z])
    #print(coordxy)


    #coordxy[0][:] = x

    return coordxy

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
    r = 0.15
    coord_size = len(coords)


    #Main Loop which iterates over all points and calculates A
    for id in range(coord_size):

        my_list = create_neighborhood(r, coords,id)

        #print('mylist', my_list)
        #print(len(my_list))
        #print(coords[my_list[id]])

        #A = constructa(coords, my_list)

        #cnd = linalg.cond(A)
        #print('Condition number of A ', cnd)

        #InvA = inverta(A)
        p = projection(coords,my_list)

        #print('P Matrix:', p)
        print('***End of loop ', id, '***')

        end = time.time()

    print('time elapsed:', end-start)
#print(A)

if __name__ == '__main__':
    print('Running read points buildA file')
    execute_main()
