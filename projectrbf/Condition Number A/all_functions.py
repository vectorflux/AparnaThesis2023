import numpy as np
import netCDF4 as nc
import scipy.linalg as scp
import numpy.linalg as linalg
from WuWendland_functions import *
import math

def create_neighborhood(r, coords, my_i):
    coord_size = len(coords)
    list = []
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

def read_data_netcdf():

    f = nc.Dataset('/home/dlaparna/thesis1_git/MyThesis2023/projectrbf/GRID_FILES/grid.nc')
    #print(f.variables.keys())


    clon = f.variables['clon']
    clat = f.variables['clat']
    lon = clon[:]
    lat = clat[:]

    ## Convert lon lat from radians to degres. Uncomment if needed
    lon[:] = lon[:]*(180/math.pi)
    lat[:] = lat[:]*(180/math.pi)


    # R = 6371229 Radius of earth in meters
    R = 1 # radius of unit sphere
    x = R * np.cos(lat) * np.cos(lon)
    y = R * np.cos(lat) * np.sin(lon)
    z = R * np.sin(lat)

    lonlat = np.column_stack([lon,lat]) #lon lat in degrees

    coordxyz = np.column_stack([x,y,z]) #cartesian coordinates

    return coordxyz

def eucl_norm(x1,x2):
    return np.linalg.norm(x1-x2)
def getneighcoords(my_list, coords):
    xyz_r = np.zeros([len(my_list),3])
    xyz_r[:] = coords[my_list[:]]

    return xyz_r

def constructA(xyz_r):

    n = len(xyz_r)  # gets the size of the matrix
    print("size of A : ", n)
    A = np.zeros((n, n))  # declares a matrix with zeros for our interpolant matrix A

        # diagonal elements


    # elements below and above diagonal
    for i in range(n):
        for k in range(n):
            #print("Entering 2nd A for loop")
            r = eucl_norm(xyz_r[i], xyz_r[k])  # sends two tuples and gets the norm back
            #print(r)
            A[i][k] = wendland0(r)
            A[k][i] = A[i][k]

    for i in range(n):
        A[i][i] = wendland0(0)

        #print("Matrix A: ", A)
        invA = inverta(A)

        return A, invA


def inverta(A):
    p, l, u = scp.lu(A)
    inv_a = np.matmul(linalg.inv(u), linalg.inv(l))
    return inv_a


def projection(xyz_r):

    #I -xx*U
    n = len(xyz_r)
    px = np.ndarray([n, 3])
    py = np.ndarray([n, 3])
    pz = np.ndarray([n, 3])



    for i in range(n):

        x = xyz_r[i,0]
        y = xyz_r[i,1]
        z = xyz_r[i,2]


        px[i, 0] = (1-(x**2))
        px[i, 1] = (0-(x*y))
        px[i, 2] = (0-(x*z))


        py[i, 0] = (0-(x*y))
        py[i, 1] = (1-(y**2))
        py[i, 2] = (0-(y*z))


        pz[i, 0] = (0-(x*z))
        pz[i, 1] = (0-(y*z))
        pz[i, 2] = (1-(z**2))

    return px, py, pz

def differentiation (inv_A,xyz_r,Xj):

    # B = [ xj.(Xj*.Xk) -xk ]*{phi'(r(Xj))/r(Xj)}
    # inv_A = inverse of A matrix calculated before
    # Xj = xyz coordinates where the neighborhood is calculated
    # xk are the list of coordinates in the neighborhood
    # xyz_r = the coordinate list of all neighborhood points
    # D = B * inv_A

    n = len(xyz_r)
    Bx = np.zeros(n)
    By = np.zeros(n)
    Bz = np.zeros(n)
    r_x = np.zeros(n)
    mult_constant = np.empty(n)

    for i in range(n):
        r_x[i] = eucl_norm(xyz_r[i],Xj)
        mult_constant[i] = wendland0_prime(r_x[i]) / r_x[i]

    for i in range(n):
        Bx[i] = (Xj[0] * (np.matmul(Xj, np.transpose(xyz_r[i])) - np.transpose(xyz_r[i, 0]))) * mult_constant[i]
        By[i] = (Xj[1] * (np.matmul(Xj, np.transpose(xyz_r[i])) - np.transpose(xyz_r[i, 1]))) * mult_constant[i]
        Bz[i] = (Xj[2] * (np.matmul(Xj, np.transpose(xyz_r[i])) - np.transpose(xyz_r[i, 2]))) * mult_constant[i]
    ## Dimensions ##

    #[1,nrj] = np.matmul( [1,nrj],[nrj,nrj] )
    Dx = np.matmul(Bx,inv_A)
    Dy = np.matmul(By,inv_A)
    Dz = np.matmul(Bz,inv_A)

    #print('Dx: ', Dx)
    #print('Dy: ', Dy)
    #print('Dz: ', Dz)

    return Dx, Dy, Dz
