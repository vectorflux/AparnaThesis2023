import numpy as np
from WuWendland_functions import *
from rbf_matrices import *



def projection(coords, list):

    #I -xx*U
    n = len(list)
    p = np.ndarray([(n*3),3])
    coordarray = []
    #myrange = n*3
    #print('size of p:',np.shape(p))
    k = 0

    for i in range(n):

        x = coords[list[i]][0]
        y = coords[list[i]][1]
        z = coords[list[i]][2]

        #print('coords:',x,y,z )
        #print('Radius:', (x**2 + y**2 + z**2))

        #check if points on unit sphere
        #if(x**2 + y**2 +z**2 == 1):
        #    print('ok')
        #else:
        #    break

        p[k+0, 0] = (1-(x**2))
        p[k+0, 1] = (0-(x*y))
        p[k+0, 2] = (0-(x*z))
        p[k+1, 0] = (0-(x*y))
        p[k+1, 1] = (1-(y**2))
        p[k+1, 2] = (0-(y*z))
        p[k+2, 0] = (0-(x*z))
        p[k+2, 1] = (0-(y*z))
        p[k+2, 2] = (1-(z**2))
        k = k+3
        coordarray.append(x)
        coordarray.append(y)
        coordarray.append(z)


    #print('coordarray:',coordarray)
    #np.reshape(coordarray,(n*3,3))
    dotproduct = np.dot(np.transpose(coordarray),p)
    #print('Dot product is ', dotproduct)

    return p

def differentiation (coords, inv_A ,Xk):

    # B = [ xj.(Xj*.Xk) -xk ]*{phi'(r(Xj))/r(Xj)}
    # index j : iterates over all the list of subdomain
    # index k : index of the node where the calculation is being done

    # D = B * inv_A
    n = len(coords)
    Bx = np.empty(n)
    By = np.empty(n)
    Bz = np.empty(n)

    ### Loops for Bx By Bz

    for j in range(n):
        r_x = eucl_norm(coords[j],Xk)
        phi_p = wendland0_prime(r_x)

        #print('xj:', coords[j,0])
        #print('Xj:', coords[j])
        #print('Xk:', Xk)

        Bx[j] = (coords[j,0]*(np.dot(coords[j],np.transpose(Xk)))) - Xk[0]*(phi_p/r_x)

    for j in range(n):
        r_x = eucl_norm(coords[j], Xk)
        phi_p = wendland0_prime(r_x)
        By[j] = (coords[j,1]*(np.dot(coords[j],np.transpose(Xk)))) - Xk[1]*(phi_p/r_x)

    for j in range(n):
        r_x = eucl_norm(coords[j], Xk)
        phi_p = wendland0_prime(r_x)
        Bz[j] = (coords[j,2]*(np.dot(coords[j],np.transpose(Xk)))) - Xk[2]*(phi_p/r_x)


    Dx = np.matmul(Bx,inv_A)
    Dy = np.matmul(By,inv_A)
    Dz = np.matmul(Bz,inv_A)

    print('Dx: ', Dx)
    print('Dy: ', Dy)
    print('Dz: ', Dz)












