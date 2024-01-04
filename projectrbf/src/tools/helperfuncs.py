import numpy as np
import math
#from initializefields import *
from testcases import *


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


def gvec2cvec_sphere(lon_fun, lat_fun, lam, th,ghost):

    
    #lat_fun = lat_fun[:,np.newaxis]
    #lon_fun = lon_fun[:,np.newaxis]


    l = len(lat_fun)

    #vlat = np.tile(lat_fun,3)
    #vlon = np.tile(lon_fun,3)

    #s2c_lat = np.zeros((l,3))
    #s2c_lon = np.zeros((l,3))

    uvw = np.zeros((l,3))

    for i in range(l):
        
            #s2c_lat[i,0] = -np.sin(lam[i]) 
            #s2c_lat[i,1] = np.cos(lam[i]) 
            #s2c_lat[i,2] = 0.0
            #s2c_lon[i,0] = -np.cos(lam[i])*np.sin(th[i]) 
            #s2c_lon[i,1] = -np.sin(lam[i])*np.sin(th[i]) 
            #s2c_lon[i,2] = np.cos(th[i])
    
            s2c_u = np.array([-np.sin(lam[i]), np.cos(lam[i]), 0.0])
            s2c_v = np.array([-np.cos(lam[i])*np.sin(th[i]), -np.sin(lam[i])*np.sin(th[i]), np.cos(th[i])])
            uvw[i,0] = s2c_u[0]*lat_fun[i] + s2c_v[0]*lon_fun[i]
            uvw[i,1] = s2c_u[1]*lat_fun[i] + s2c_v[1]*lon_fun[i]
            uvw[i,2] = s2c_u[2]*lat_fun[i] + s2c_v[2]*lon_fun[i]

    
    #print("Shape of vlat:", np.shape(vlat))
    #print("Shape of s2c:", np.shape(s2c_lon))
    #print(" uvw:", uvw)
    return uvw


def cvec2gvec_sphere(uvwh, lon, lat):    

    l = len(lon)
    
    #print("Lengths of uvwh and lon, lat: ", len(uvwh),len(lon),len(lat) )
    p_gu = np.zeros(l)
    p_gv = np.zeros(l)
    p_norm = np.zeros(l)

    for i in range(l):
        z_sln = np.sin(lon[i])
        z_cln = np.cos(lon[i])
        z_slt = np.sin(lat[i])
        z_clt = np.cos(lat[i])

        p_gu[i] = (z_cln*uvwh[i,1]) - (z_sln*uvwh[i,0])
        tmp = (z_cln*uvwh[i,0]) + (z_sln*uvwh[i,1])
        tmp = (z_slt*tmp)
        p_gv[i] = (z_clt*uvwh[i,2]) - tmp
        
        #v_lon, v_lat, h = test2_fun(lon[i],lat[i])

        #if p_gu[i]==v_lat and p_gv[i]==v_lon:
            #print("yay")
        #else:
            #print("***Nay***")

        #v_lon[n], v_lat[n], h[n] = test6_fun(lam[n],th[n])
        
        #p_norm[i] = np.sqrt(p_gu[i]**2 + p_gv[i]**2) 

        #p_gu[i] = p_gu[i]/p_norm[i]
        #p_gv[i] = p_gv[i]/p_norm[i]

    #for i in range(l):
        
        #z_sln = np.sin(lon[i])
        #z_cln = np.cos(lon[i])
        #z_slt = np.sin(lat[i])
        #z_clt = np.cos(lat[i])

        #c2s_u = np.array([-z_sln, -z_slt*z_cln])
        #c2s_v = np.array([z_cln, -z_slt*z_sln])
        #c2s_w = np.array([0.0, z_clt])

        #p_gu[i] = (c2s_u[0]*uvwh[i,0]) + (c2s_v[0]*uvwh[i,1]) + (c2s_w[0]*uvwh[i,2])
        #p_gv[i] = (c2s_u[1]*uvwh[i,0]) + (c2s_v[1]*uvwh[i,1]) + (c2s_w[1]*uvwh[i,2])

    #print("GeoVelocity:", geovel)
    return p_gu, p_gv

def validate_halo_exchange(uvwh,xyz, n_p,ghost):
    vt = 1
    for n in range(n_p):
        if ghost[n]:
            if (uvwh[n,3] != 4*(xyz[n,0]**2 +xyz[n,1]**2 +xyz[n,2]**2)):
                print("Halo exchange failed")
                vt = 0

    return vt

