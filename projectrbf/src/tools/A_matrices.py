############
# INPUT : An 2D array of coordinates which form the neighborhood cluster for the node at which we are calculating
# OUTPUT : Returns A (interpolation matrix) and invA (inverse of interpolation matrix)
# WORKING : Can change function call statement in line #27 and line #37 for a different Wendland function. Check the file WuWendland_functions.py 

#############

import numpy as np
import scipy.linalg as scp
import numpy.linalg as linalg
from WuWendland_functions import *
from helperfuncs import *


def constructA(xyz_r,normalization_factor):

        n = len(xyz_r) #gets the size of the matrix
        A = np.zeros((n, n))  #declares a matrix with zeros for the construction of our interpolant matrix A


        # diagonal elements
        for i in range(n):
            A[i][i] = wendland1(0) #can change as per need

        #elements below and above diagonal
        for i in range(n):
            for k in range(n):
                if k < i :
                    r = eucl_norm(xyz_r[i],xyz_r[k]) #sends two tuples and gets the norm back
                    r=r/normalization_factor  #scaling to match the Wendland functions
#                    A[i][k] = wendland1(r)*normalization_factor #can change as per need
                    A[i][k] = wendland1(r) #can change as per need
                    A[k][i] = A[i][k]

        invA = inverta(A)

        return A, invA



def inverta(A):

    pl,u = scp.lu(A, permute_l = True)
    u_inv = linalg.inv(u)
    pl_inv = linalg.inv(pl)

    inv_a = np.matmul(u_inv,pl_inv)
    
    return inv_a






