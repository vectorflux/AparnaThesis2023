
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
    #r_x = np.zeros(n)
    #mult_constant = np.empty(n)
    
    for i in range(n):
        #print(xyz_r[i,:])
        r_x = eucl_norm(xyz_r[i,:],Xj)/0.065 #scaling to match wendland functions
        mult_constant = wendland1_prime(r_x) / r_x

 
    #[1,nrj] = ([1,1]*(np.matmul([1,3],[3,nrj]) - [1,nrj]) * [1,nrj]
        Bx[i] = ((xyz_r[i,0]*(np.dot(Xj,(xyz_r[i,:])))) - Xj[0])*mult_constant
        By[i] = ((xyz_r[i,1]*(np.dot(Xj,(xyz_r[i,:])))) - Xj[1])*mult_constant
        Bz[i] = ((xyz_r[i,2]*(np.dot(Xj,(xyz_r[i,:])))) - Xj[2])*mult_constant
        
        Bxyz = np.row_stack([Bx[i],By[i],Bz[i]])
        product = np.dot(xyz_r[i,:],Bxyz)
        
        #print("dotproduct of B:", product)

    #[1,nrj] = np.matmul( [1,nrj],[nrj,nrj] )

    Bx = Bx[np.newaxis,:]
    By = By[np.newaxis,:]
    Bz = Bz[np.newaxis,:]
    
    Bxyz = np.row_stack([Bx,By,Bz])
    #print("Bx:", Bxyz)
    
    Dx = np.matmul(Bx,inv_A)
    Dy = np.matmul(By,inv_A)
    Dz = np.matmul(Bz,inv_A)

    Dxyz = np.row_stack([Dx,Dy,Dz])
    
    #print(np.shape(Xj), np.shape(Dxyz))
    #dotproduct = np.dot(Xj,Dxyz)
    #for i in range(len(xyz_r)):
        #product = np.dot(xyz_r[i,:],Dxyz[:,i])
        #print( " dotproduct  of dxyz:", product)

    #print('Dx: ', Dx)
    #print('Dy: ', Dy)
    #print('Dz: ', Dz)

    #print("Shapes of Dx:", np.shape(Dx))

    return Dx, Dy, Dz















