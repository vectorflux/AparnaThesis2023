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
            A[i][i] = wendland0(0)

        #elements below and above diagonal
        for i in range(n):
            for k in range(n):
                if k < i :
                    r = eucl_norm(xyz_r[i],xyz_r[k]) #sends two tuples and gets the norm back
                    #print(r)
                    A[i][k] = wendland0(r)
                    A[k][i] = A[i][k]

        invA = inverta(A)

        np.reshape(A, (1,nrj))
        np.reshape(invA, (1,nrj))


        return A, invA



def inverta(A):

    p,l,u = scp.lu(A)
    inv_a = np.matmul(linalg.inv(u),linalg.inv(l))
    return inv_a


#def storeallA():





