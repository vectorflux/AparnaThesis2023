#calls read netcdf file
# relays data to atlasfuncinterface
# gets indices back

from read_netcdf import *
from atlas_func_interface import *

lonlat = read_data_netcdf()

print(lonlat)

myradius = 0.065

indices = find_neighbors(myradius,lonlat)

print('Indices: ', indices)