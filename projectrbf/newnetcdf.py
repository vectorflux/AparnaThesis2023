import netCDF4 as nc


f = nc.Dataset('/users/ldevulap/MyThesis2023/projectrbf/GRID_FILES/grid.nc')
print(f.variables.keys())
