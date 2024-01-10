
import numpy as np
import scipy.linalg as scp
import numpy.linalg as linalg
from WuWendland_functions import *
from helperfuncs import *


def constructA(xyz_r):

        n = len(xyz_r) #gets the size of the matrix
        A = np.zeros((n, n))  #declares a matrix with zeros for our interpolant matrix A


        # diagonal elements
        for i in range(n):
            A[i][i] = wendland1(0)

        #elements below and above diagonal
        for i in range(n):
            for k in range(n):
                if k < i :
                    r = eucl_norm(xyz_r[i],xyz_r[k]) #sends two tuples and gets the norm back
                    #print(r)
                    r=r/0.065 #scaling to match the Wendland functions
                    #print("r:", r)
                    A[i][k] = wendland1(r)
                    A[k][i] = A[i][k]

        invA = inverta(A)
        
        #det = linalg.det(np.dot(A,invA))
        #if not np.isclose(det,1.0):
            #print(det)

        #nrj = n*n
        #A = np.reshape(A, (1,nrj))
        #np.reshape(invA, (1,nrj))


        return A, invA



def inverta(A):

    pl,u = scp.lu(A, permute_l = True)
    #l = np.dot(p,l)

    u_inv = linalg.inv(u)


    pl_inv = linalg.inv(pl)
    
    #print(np.allclose(A, pl@u))

    inv_a = np.matmul(u_inv,pl_inv)
    return inv_a


#def storeallA():





