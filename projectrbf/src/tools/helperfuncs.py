import numpy as np
import math
#from initializefields import *
from testcases import *




def getpxyz(xyz):

    px = np.zeros(3)
    py = np.zeros(3)
    pz = np.zeros(3)

    x = xyz[0]
    y = xyz[1]
    z = xyz[2]


        #if i==1:
            #print("xyz", x,y,z)

    px[0] = (1 - (x ** 2))
    px[1] = (0 - (x * y))
    px[2] = (0 - (x * z))

    py[0] = (0 - (x * y))
    py[1] = (1 - (y ** 2))
    py[2] = (0 - (y * z))

    pz[0] = (0 - (x * z))
    pz[1] = (0 - (y * z))
    pz[2] = (1 - (z ** 2))

    px = px[np.newaxis,:] #row matrix
    py = py[np.newaxis,:]
    pz = pz[np.newaxis,:]

    return px, py, pz
def get_uvwh_r(uvwh,nearest):
    uvwh_r = np.zeros([len(nearest),4])
    for n in range(len(nearest)):
        uvwh_r[n,0] = uvwh[int(nearest[n]),0]
        uvwh_r[n,1] = uvwh[int(nearest[n]),1]
        uvwh_r[n,2] = uvwh[int(nearest[n]),2]
        uvwh_r[n,3] = uvwh[int(nearest[n]),3]

        #print(np.where((np.sum(uvwh_r, axis=1))==0)[0])


    return uvwh_r



def getcartesian(lonlat):

    #lat = lonlat[:,1]
    #lon = lonlat[:,0]
    

    ## Convert lon lat from degrees to radians. Uncomment if needed
    lon = np.radians(lonlat[:,0])
    lat = np.radians(lonlat[:,1])

    # R = 6371229 Radius of earth in meters
    R = 1.0  # radius of unit sphere
    x = R * np.cos(lat) * np.cos(lon)
    y = R * np.cos(lat) * np.sin(lon)
    z = R * np.sin(lat)

    xyz = np.column_stack([x,y,z])

    return xyz


def eucl_norm(x1,x2):
    return np.linalg.norm(x1-x2)

def getneighcoords(my_list, coords):
    xyz_r = np.zeros((len(my_list),3))
    #print("My list:", my_list) 
    for i in range(len(my_list)):
        k = int(my_list[i])
        xyz_r[i,:] = coords[k]

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
        
        #v_lon, v_lat, h,f = test2_fun(lon[i],lat[i])
        
        #diff_u = np.abs(v_lat - p_gu[i])
        #diff_v = np.abs(v_lon - p_gv[i])

        #p_gu[i] = v_lat
        #p_gv[i] = v_lon

        #if diff_u <= 10e-12:
            #print("yay for u")
        #else:print("wrong for u")
        
        #if diff_v<= 10e-12:
            #print("yay for v")

        #else:
            #print("wrong for v")

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

