import numpy as np
from A_matrices import *
from WuWendland_functions import *



def differentiation (inv_A,xyz_r,Xj,normalization_factor):


    # B = [ xk.(Xj*.Xk) -xj ]*{phi'(r(Xj))/r(Xj)}
    # inv_A = inverse of A matrix calculated before
    # Xj = xyz coordinates where the neighborhood is calculated
    # xyz_r = the coordinate list of all neighborhood points
    # D = B * inv_A

    n = len(xyz_r)
    Bx = np.zeros(n)
    By = np.zeros(n)
    Bz = np.zeros(n)
    
    for i in range(n):


        #r_x = np.linalg.norm(xyz_r[i,:]-Xj)/normalization_factor #scaling to match wendland functions
        #mult_constant = wendland_3_1_prime(r_x) / r_x / normalization_factor # derivative

        r_x = eucl_norm(xyz_r[i,:],Xj)/normalization_factor #scaling to match wendland functions
        mult_constant = wendland1_prime(r_x) / r_x /normalization_factor


 
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

def gradient(c,xyz_r,Xj,normalization_factor):

    # Determine the unprojected 3D gradient

    # G = [ xk -xj ]*{phi'(r(Xj))/r(Xj)}
    # Xj = xyz coordinates where the neighborhood is calculated
    # xyz_r = the coordinate list of all neighborhood points
    # normalization_factor = the actual cutoff radius

    n = len(xyz_r)
    G = np.zeros(3)
    for i in range(n):
        r_x = np.linalg.norm(xyz_r[i,:]-Xj)/normalization_factor #scaling to match wendland functions
        mult_constant = wendland_3_1_prime(r_x) / r_x / normalization_factor  # derivative
        G += c[i]*(xyz_r[i,:] - Xj)*mult_constant
        
    return G

def proj_gradient(c,xyz_r,Xj,normalization_factor):

    # Determine the unprojected 3D gradient

    # G = [ xk -xj ]*{phi'(r(Xj))/r(Xj)}
    # Xj = xyz coordinates where the neighborhood is calculated
    # xyz_r = the coordinate list of all neighborhood points
    # normalization_factor = the actual cutoff radius

    n = len(xyz_r)
    B = np.zeros(3)
    for i in range(n):
        r_x = np.linalg.norm(xyz_r[i,:]-Xj)/normalization_factor #scaling to match wendland functions
        mult_constant = wendland_3_1_prime(r_x) / r_x / normalization_factor # derivative
        B[0] += c[i]*((xyz_r[i,0]*(np.dot(Xj,(xyz_r[i,:])))) - Xj[0])*mult_constant
        B[1] += c[i]*((xyz_r[i,1]*(np.dot(Xj,(xyz_r[i,:])))) - Xj[1])*mult_constant
        B[2] += c[i]*((xyz_r[i,2]*(np.dot(Xj,(xyz_r[i,:])))) - Xj[2])*mult_constant

    return B
