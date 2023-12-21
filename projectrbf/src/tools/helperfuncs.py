import numpy as np
import math


def getcartesian(lonlat):

    #lat = lonlat[:,1]
    #lon = lonlat[:,0]
    

    ## Convert lon lat from degrees to radians. Uncomment if needed
    lon = np.radians(lonlat[:,0])
    lat = np.radians(lonlat[:,1])

    # R = 6371229 Radius of earth in meters
    R = 1.  # radius of unit sphere
    x = R * np.cos(lat) * np.cos(lon)
    y = R * np.cos(lat) * np.sin(lon)
    z = R * np.sin(lat)

    xyz = np.column_stack([x,y,z])

    return xyz


def eucl_norm(x1,x2):
    return np.linalg.norm(x1-x2)

def getneighcoords(my_list, coords):
    xyz_r = np.zeros([len(my_list),3])
    xyz_r[:] = coords[my_list[:]]

    return xyz_r


def cartesian_fromlonlat(lat_fun, lon_fun, lam, th):

    l = len(lat_fun)
    
    #uvw = np.zeros([l,3])
    
    lat_fun = lat_fun[:,np.newaxis]
    lon_fun = lon_fun[:,np.newaxis]

    vlat = np.tile(lat_fun,(1,3))
    vlon = np.tile(lon_fun,(1,3))

    s2c_lat = np.empty([l,3])
    s2c_lon = np.empty([l,3])

    for i in range(l):
        s2c_lat[i,0] = -np.sin(lam[i]) 
        s2c_lat[i,1] = np.cos(lam[i]) 
        s2c_lat[i,2] = 0.0


        s2c_lon[i,0] = -np.cos(lam[i])*np.sin(th[i]) 
        s2c_lon[i,0] = -np.sin(lam[i])*np.sin(th[i]) 
        s2c_lon[i,0] = np.cos(th[i])
    
    uvw = np.multiply(vlat,s2c_lat) + np.multiply(vlon,s2c_lon)
    
    return uvw 

def validate_halo_exchange(uvwh,xyz, n_p,ghost):
    vt = 1
    for n in range(n_p):
        if ghost[n]:
            if (uvwh[n,3] != 4*(xyz[n,0]**2 +xyz[n,1]**2 +xyz[n,2]**2)):
                print("Halo exchange failed")
                vt = 0

    return vt

