import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import math
from helperfuncs import *


def plot_global(uvwh, lonlat):

    lam = np.radians(lonlat[:,0])
    th = np.radians(lonlat[:,1])

    
    N = len(lam)
    xyz = getcartesian(lonlat)
    
    uvw = np.zeros((N,3))
    h = np.zeros(N)
    v_lon = np.zeros(N)
    v_lat = np.zeros(N)
    
    #print("length of xyz: ", len(xyz))
    x = xyz[:,0]
    y = xyz[:,1]
    z = xyz[:,2]
    
    u = uvwh[:,0]
    v = uvwh[:,1]
    w = uvwh[:,2]
    h = uvwh[:,3]

    v_zonal, v_meri = cvec2gvec_sphere(uvwh,lam,th)
    
    # Set up a 3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    scatter = ax.scatter(x,y,z, c=h, cmap='jet')
    fig.colorbar(scatter, ax=ax, label='Values')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('Geopotential height' )  
    plt.savefig('height.png')
    #plt.show()
    plt.close()

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    scatter2 = ax.scatter(x,y,z, c=v_zonal, cmap='jet')
    fig.colorbar(scatter2, ax=ax, label='Values')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('Zonal velocity' )
    plt.savefig('zonal.png')
    #plt.show()
    plt.close()
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    scatter3 = ax.scatter(x,y,z, c=v_meri, cmap='jet')
    fig.colorbar(scatter3, ax=ax, label='Values')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('Meridional Velocity' )
    plt.savefig('meri.png')
    #plt.show()
    plt.close()



def plot_global2(field,lonlat):

    lam = np.radians(lonlat[:,0])
    th = np.radians(lonlat[:,1])


    N = len(lam)
    xyz = getcartesian(lonlat)
 
    #print("length of xyz: ", len(xyz))
    x = xyz[:,0]
    y = xyz[:,1]
    z = xyz[:,2]

    c_zonal, c_meri = cvec2gvec_sphere(field,lam,th)
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    scatter = ax.scatter(x,y,z, c=c_zonal, cmap='jet')
    fig.colorbar(scatter, ax=ax, label='Values')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('Termc - Zonal height update' )
    plt.savefig('termc_zonal.png')
    #plt.show()
    plt.close()

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    scatter2 = ax.scatter(x,y,z, c=c_meri, cmap='jet')
    fig.colorbar(scatter, ax=ax, label='Values')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('Termc - Meridional height update' )
    plt.savefig('termc_meri.png')
    #plt.show()
    plt.close()




