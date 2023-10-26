
import netCDF4 as nc
import numpy as np
import math

def read_data_netcdf():

    f = nc.Dataset('/users/ldevulap/GRID_FILES/grid.nc')
    #print(f.variables.keys())


    clon = f.variables['clon']
    clat = f.variables['clat']
    lon = clon[:]
    lat = clat[:]

    ## Convert lon lat from radians to degres. Uncomment if needed
    lon[:] = lon[:]*(180/math.pi)
    lat[:] = lat[:]*(180/math.pi)


    # R = 6371229 Radius of earth in meters
    R = 1 # radius of unit sphere
    x = R * np.cos(lat) * np.cos(lon)
    y = R * np.cos(lat) * np.sin(lon)
    z = R * np.sin(lat)


    lonlat = np.column_stack([lon,lat])


    return lonlat
