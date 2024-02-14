
import numpy as np
#from operator_matrices import *
import math
from helperfuncs import *
#from testcases import *



def construct_rhsd_vector(nrj_size_list, allnearest, normalization_factor, uvwh, f, xyz, allD, ghost):
    n_p = len(uvwh)
    rhsd = np.zeros(3)
    it = 0
    termA = np.zeros(3)
    termB = np.zeros(3)
    termC = np.zeros(3)

    R  = np.zeros((n_p,4))

    g = 1.5454e-6
    l = 0

    # Calculating RHS values for each of the internal points
    for i in range(n_p):
        if not ghost[i]:

            k = int(nrj_size_list[l]) #Check the requirement of int typecasting here
            l= l+1
            Dnx, Dny, Dnz, nearest = get_Dnxyz_from_allD(allD,allnearest, it,k)

            it = it + k

            xyz_r = getneighcoords(nearest, xyz)
            uvwh_r = get_uvwh_r(uvwh, nearest)

            u = uvwh_r[:, 0]
            v = uvwh_r[:, 1]
            w = uvwh_r[:, 2]
            h = uvwh_r[:, 3]

            termA[0] = (uvwh[i,0]*(np.dot(np.transpose(Dnx),u)) +
                 uvwh[i,1]*(np.dot(np.transpose(Dny),u)) +
                 uvwh[i,2]*(np.dot(np.transpose(Dnz),u)))

            termA[1] = (uvwh[i, 0] * (np.dot(np.transpose(Dnx), v)) +
                       uvwh[i, 1] * (np.dot(np.transpose(Dny), v)) +
                       uvwh[i, 2] * (np.dot(np.transpose(Dnz), v)))

            termA[2] = (uvwh[i, 0] * (np.dot(np.transpose(Dnx), w)) +
                    uvwh[i, 1] * (np.dot(np.transpose(Dny), w)) +
                    uvwh[i, 2] * (np.dot(np.transpose(Dnz), w)))

            #f. [3,1]
            termB[0] = f[i]*((xyz[i,1]*uvwh[i,2]) - (xyz[i,2]*uvwh[i,1]))
            termB[1] = f[i]*((xyz[i,2]*uvwh[i,0]) - (xyz[i,0]*uvwh[i,2]))
            termB[2] = f[i]*((xyz[i,0]*uvwh[i,1]) - (xyz[i,1]*uvwh[i,0]))

            #g.[Dnx Dny Dnz]. h
            termCx= np.transpose(Dnx)
            termCy= np.transpose(Dny)
            termCz= np.transpose(Dnz)

            term_c = np.row_stack([termCx,termCy,termCz])
            termC = g*(np.dot(term_c,h))

            rhsd[0] = termA[0] + (termB[0]) + (termC[0])
            rhsd[1] = termA[1] + (termB[1]) + (termC[1])
            rhsd[2] = termA[2] + (termB[2]) + (termC[2])

            #Get px, py, pz at the node where we are calculating
            px, py, pz = getpxyz(xyz[i,:])


            R[i,0] = -np.dot(px, rhsd)
            R[i,1] = -np.dot(py, rhsd)
            R[i,2] = -np.dot(pz, rhsd)

            R[i,3] = (uvwh[i, 0]*(np.dot(np.transpose(Dnx),h)) + uvwh[i, 1]*(np.dot(np.transpose(Dny),h)) + uvwh[i, 2]*(np.dot(np.transpose(Dnz),h)) + uvwh[i, 3]*(np.dot(np.transpose(Dnx),u) + np.dot(np.transpose(Dny),v)+ np.dot(np.transpose(Dnz),w)))

    return R

