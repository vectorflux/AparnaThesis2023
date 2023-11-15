
import numpy as np
from A_matrices import *
from WuWendland_functions import *



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

    ## Dimensions ##


    #[1,nrj] = ([1,1]*(np.matmul([1,3],[3,nrj]) - [1,nrj]) * [1,nrj]
    Bx[:] = (Xj[0]*(np.matmul(Xj, np.transpose(xyz_r[:])) - np.transpose(xyz_r[:,0])))*mult_constant[:]
    By[:] = (Xj[1]*(np.matmul(Xj, np.transpose(xyz_r[:])) - np.transpose(xyz_r[:,1])))*mult_constant[:]
    Bz[:] = (Xj[2]*(np.matmul(Xj, np.transpose(xyz_r[:])) - np.transpose(xyz_r[:,2])))*mult_constant[:]

    ## OLD 2 ## Bx[:] = (xyz_r[:,0]*(np.dot(xyz_r[j],np.transpose(Xj))) - Xj[0])*mult_constant[:]

    ## OLD ##By[:] = (coords[j,1]*(np.dot(coords[j],np.transpose(Xk)))) - Xk[1]*(phi_p/r_x)

    #[1,nrj] = np.matmul( [1,nrj],[nrj,nrj] )
    Dx = np.matmul(Bx,inv_A)
    Dy = np.matmul(By,inv_A)
    Dz = np.matmul(Bz,inv_A)

    #print('Dx: ', Dx)
    #print('Dy: ', Dy)
    #print('Dz: ', Dz)



    return Dx, Dy, Dz















