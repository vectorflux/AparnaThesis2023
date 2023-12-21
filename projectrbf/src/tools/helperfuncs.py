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
