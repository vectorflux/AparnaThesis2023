import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import math
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from helperfuncs import *
from read_netcdf_file import *
from initializefields import *


def plot_global(uvwh, lonlat):

    #lonlat = read_data_netcdf()
    
    lon =lonlat[:,0]
    lat =lonlat[:,1]

    lam = np.radians(lon[:])
    th = np.radians(lat[:])

    #th = math.acos(z/(math.sqrt(x**2 +y**2 + z**2)))
    ## Convert lon lat from radians to degres. Uncomment if needed
    #lon[:] = lon[:]*(180/math.pi)
    #lat[:] = lat[:]*(180/math.pi)

    N = len(lam)
    xyz = getcartesian(lonlat)
    
    uvw = np.zeros([N,3])
    h = np.zeros(N)
    v_lon = np.zeros(N)
    v_lat = np.zeros(N)
    
    #print("length of xyz: ", len(xyz))
    x = xyz[:,0]
    y = xyz[:,1]
    z = xyz[:,2]
    
    
    for n in range(N):

        #v_lon[n], v_lat[n], h[n] = test2_fun(lam[n],th[n])
        v_lon[n], v_lat[n], h[n] = test6_fun(lam[n],th[n])


    #zonal_fun = uu0 * (np.cos(th) * np.cos(angle) + np.sin(th) * np.cos(lam) * np.sin(angle))
    #meri_fun = -uu0 * np.sin(angle) * np.sin(lam)
    #h0_fun =  h0 - 1 / g * ( Omega * uu0 + uu0 ** 2 / 2) * (np.sin(th) * np.cos(angle) - np.cos(lam) * np.cos(th) * np.sin(angle)) ** 2

    #print("Shape of LON-LAT:", len(lon),len(lat),"\n")
    #u = uvwh[:,0]
    #v = uvwh[:,1]
    #w = uvwh[:,2]
    #h = uvwh[:,3]

    # Set up a 3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    scatter = ax.scatter(x,y,z, c=v_lon, cmap='plasma')
    #plt.scatter(lam, th, cmap='plasma')
    #ax.quiver(x, y, z, u, v, w, length=0.1, normalize=True, color='blue', arrow_length_ratio=0.3)
    fig.colorbar(scatter, ax=ax, label='data')

    #ax = fig.add_subplot(1,1,1,projection = ccrs.Orthographic(central_latitude=0, central_longitude=0))# Plot the surface with color representing the magnitude of velocity
    #scatter = ax.scatter(lon,lat,c=z,cmap='RdBu_r',s =1,vmin=1e-6, vmax=9e-6,transform=ccrs.PlateCarree())
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    #plt.xlabel('X')
    #plt.ylabel('Y')
    #ax.coastlines (resolution='110m', color='k')
    #ax.set_global()
    ax.set_title('Meridinial Velocity on Spherical Surface' )
    #ax.set_facecolor('gray')
    #fig.colorbar(scatter, ax=ax)
    #plt.subplots_adjust(hspace=0.5)
    #plt.show()
    
    # Add a color bar
    #fig.colorbar(scatter, ax=ax, shrink=0.5, aspect=10)# Set plot labels
    #ax.set_xlabel('X-axis')
    #ax.set_ylabel('Y-axis')
    #ax.set_zlabel('Z-axis')
    #ax.set_title('3D Velocity Field on Spherical Surface')
    
    plt.savefig('velocity3d.png')
    #plt.show()



def plot_global_contour(uvwh):

    
    lonlat = read_data_netcdf()

    lon =lonlat[:,0]
    lat =lonlat[:,1]

    z = np.sqrt(np.square(uvwh[:,0]) + np.square(uvwh[:,1]) + np.square(uvwh[:,2]))
    #z = uvwh[:,0]
    #z = uvwh[:,1]
    #z = uvwh[:,2]
    
    
    fig, ax = plt.subplots(subplot_kw={'projection': ccrs.Orthographic(central_longitude=0, central_latitude=0)})# Create a contour plot
    
    
    #contour = ax.contour(lon, lat,z, levels=10, cmap='rainbow', transform=ccrs.PlateCarree())# Add coastlines for reference
    
    # Create a quiver plot
    quiver = ax.quiver(lons, lats, u, v, transform=ccrs.PlateCarree(), scale=30)
    
    ax.coastlines()# Set plot labels and title
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title('Quiver Plot with Orthographic Projection')# Add gridlines with labels
    ax.gridlines(draw_labels=True, linestyle='--', linewidth=0.5, color='black')
    ax.xaxis.set_major_formatter(LongitudeFormatter())
    ax.yaxis.set_major_formatter(LatitudeFormatter())# Add a colorbar
    #fig.colorbar(contour, ax=ax, label='Contour Levels')
    

    plt.savefig('velocity3d.png')

