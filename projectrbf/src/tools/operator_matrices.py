import numpy as np
from A_matrices import *
from WuWendland_functions import *


def differentiation (inv_A,xyz_r,Xj,norm_factor):

    # B = [ xk.(Xj*.Xk) -xj ]*{phi'(r(Xj))/r(Xj)}
    # inv_A = inverse of A matrix calculated before
    # Xj = xyz coordinates where the neighborhood is calculated
    # xk are the list of coordinates in the neighborhood
    # xyz_r = the coordinate list of all neighborhood points
    # D = B * inv_A

    n = len(xyz_r)
    Bx = np.zeros(n)
    By = np.zeros(n)
    Bz = np.zeros(n)
    
    for i in range(n):

        r_x = eucl_norm(xyz_r[i,:],Xj)/norm_factor #scaling to match wendland functions
        mult_constant = wendland1_prime(r_x) / r_x /norm_factor

 
        #[1,nrj] = ([1,1]*(np.matmul([1,3],[3,nrj]) - [1,nrj]) * [1,nrj]
        Bx[i] = ((xyz_r[i,0]*(np.dot(Xj,(xyz_r[i,:])))) - Xj[0])*mult_constant
        By[i] = ((xyz_r[i,1]*(np.dot(Xj,(xyz_r[i,:])))) - Xj[1])*mult_constant
        Bz[i] = ((xyz_r[i,2]*(np.dot(Xj,(xyz_r[i,:])))) - Xj[2])*mult_constant

    #[1,nrj] = np.matmul( [1,nrj],[nrj,nrj] )
    Bx = Bx[np.newaxis,:]
    By = By[np.newaxis,:]
    Bz = Bz[np.newaxis,:]
    
    #Bxyz = np.row_stack([Bx,By,Bz])
    
    Dx = np.matmul(Bx,inv_A)
    Dy = np.matmul(By,inv_A)
    Dz = np.matmul(Bz,inv_A)

    #Dxyz = np.row_stack([Dx,Dy,Dz])
    
    return Dx, Dy, Dz















