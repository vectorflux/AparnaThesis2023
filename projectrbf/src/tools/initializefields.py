import numpy as np
import math
from helperfuncs import *
from testcases import *


def set_initial_conditions(xyz, n_p, ghost,lonlat):
    
    N = n_p #internal+halo , 5 field columns
    f = np.zeros(N)

    fields = np.zeros((N,5))
    lam = np.radians(lonlat[:,0])
    th = np.radians(lonlat[:,1])

    x = xyz[:,0]
    y = xyz[:,1]
    z = xyz[:,2]
    
    v_lon = np.zeros(N)
    v_lat = np.zeros(N)
    h = np.zeros(N)

    for n in range(N):
        if not ghost[n]: #calculating all internal values     
            v_lon[n], v_lat[n], h[n], f[n] = test0_fun(lam[n],th[n]) #can change to test2_fun() , test6_fun()

        else:
            v_lon[n] = v_lat[n] = h[n] = 0.0
            f[n] = 0.0
    
    uvw = gvec2cvec_sphere(v_lon, v_lat, lam, th,ghost)
            
    for n in range(n_p):
        if not ghost[n]:
            fields[n,0] = uvw[n,0]
            fields[n,1] = uvw[n,1]
            fields[n,2] = uvw[n,2]
            fields[n,3] = h[n]

        else:
            fields[n,0] = 0.
            fields[n,1] = 0.
            fields[n,2] = 0.
            fields[n,3] = 0.


    return fields, f
