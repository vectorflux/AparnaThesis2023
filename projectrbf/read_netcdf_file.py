# Master Thesis Project
# This python script reads data from a given NetCDF File
# Latitudes and Longitudes are given in Radians.
# Coordinates are extracted from Lon,Lat Cartesian form


# Authors: Aparna Devulapalli



import netCDF4 as nc
import numpy as np
import math

def read_data_netcdf():

    f = nc.Dataset('/home/dlaparna/thesis/projectrbf/GRID_FILES/grid.nc')
    #print(f.variables.keys())


    clon = f.variables['clon']
    clat = f.variables['clat']
    lon = clon[:]
    lat = clat[:]

    ## Convert lon lat from radians to degres. Uncomment if needed
    #lon[:] = lon[:]*(180/math.pi)
    #lat[:] = lat[:]*(180/math.pi)

    cx = f.variables['cartesian_x_vertices']
    cy = f.variables['cartesian_y_vertices']
    cz = f.variables['cartesian_z_vertices']

    # R = 6371229 Radius of earth in meters
    R = 1 # radius of unit sphere
    x = R * np.cos(lat) * np.cos(lon)
    y = R * np.cos(lat) * np.sin(lon)
    z = R * np.sin(lat)


    coordxy = np.column_stack([x,y,z])

    return coordxy


