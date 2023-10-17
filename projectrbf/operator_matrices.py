import numpy as np



def projection(coords, list):



    #I -xx*U
    n = len(list)
    p = np.ndarray([(n*3),3])
    #myrange = n*3
    print('size of p:',np.shape(p))
    k = 0

    for i in range(n):

        x = coords[list[i]][0]
        y = coords[list[i]][1]
        z = coords[list[i]][2]

        #check if points on unit sphere
        #if(x**2 + y**2 +z**2 != 1):




        p[k+0, 0] = (1-(x**2))
        p[k+0, 1] = (0-(x*y))
        p[k+0, 2] = (0-(x*z))
        p[k+1, 0] = (0-(x*y))
        p[k+1, 1] = (1-(y**2))
        p[k+1, 2] = (0-(y*z))
        p[k+2, 0] = (0-(x*z))
        p[k+2, 1] = (0-(y*z))
        p[k+2, 2] = (1-(z**2))

        dotproduct = np.dot(p,coords[list[i]])
        print('Cross product is ', crossproduct)

        k = k+3
    return p

def gradient():

    dx =
    dy =
    dz =

    np.gradient(f,axis=)

    #take a vector and returns the gradient of it in x y z direction











