
import scipy.integrate as integrate
import numpy as np
import scipy.linalg as scp
import numpy.linalg as linalg
from WuWendland_functions import *


def getsubdomain(my_list, coords):
    sub_coords = np.zeros([len(my_list),3])
    sub_coords[:] = coords[my_list[:]]

    return sub_coords
def eucl_norm(x1,x2):
    return np.linalg.norm(x1-x2)

def constructa(neighbors,xyz):

    #input coordinates(xyz)

    #iterates over each row of allindices
    fullsize = len(allindices)
    all_invA = np.zeros(np.shape(allindices))



    for id in range(fullsize):
        coords = np.trim_zeros(allindices[id],'b')
        n = len(coords) #gets the size of the matrix
        A = np.zeros((n, n))  #declares a matrix with zeros for our interpolant matrix A

        # diagonal elements
        for i in range(n):
            A[i][i] = wendland0(0)

        #elements below and above diagonal
        for i in range(n):
            for k in range(n):
                if k < i :
                    r = eucl_norm(coords[i],coords[k]) #sends two tuples and gets the norm back
                    #print(r)
                    A[i][k] = wendland0(r)
                    A[k][i] = A[i][k]

        invA = inverta(A)
        all_invA = np.vstack(A)

    return all_invA

def inverta(A):

    p,l,u = scp.lu(A)
    inv_a = np.matmul(linalg.inv(u),linalg.inv(l))
    return inv_a

if __name__ == '__main__':
    #rbfmain()

    print('Running RBF Matrices file')