#### Just for one node at 1717 #####


#!/usr/bin/env python3
import atlas4py as atlas
import numpy as np
import time

from operator_matrices import *
from atlas_func_interface import *
from read_netcdf_file import *
from A_matrices import *
from initializefields import *
from visualization import *
from rk4 import *
from construct_rhsd import *


atlas.initialize() #initializes atlas and MPI

myradius = 0.065 
lonlat = read_data_netcdf()
grid = atlas.UnstructuredGrid(lonlat[:, 0], lonlat[:, 1]) #Create Unstructured Grid

#Create functionspace + partitioning #levels in the pointcloud implementation
functionspace = atlas.functionspace.PointCloud(grid, halo_radius=myradius * 2, geometry="UnitSphere")

###Access fields as numpy arrays, if needed
lonlat = atlas.make_view(functionspace.lonlat) # longitude latitude in degrees
ghost = atlas.make_view(functionspace.ghost) # ghost: 1 for ghost, 0 for owned
partition = atlas.make_view(functionspace.partition) #partition owning point (0-based)
remote_index = atlas.make_view(functionspace.remote_index) # local index on partition owning point (careful, 1-based when atlas_HAVE_FORTRAN)
global_index = atlas.make_view(functionspace.global_index) # global index across partitions (always 1-based)
xyz = getcartesian(lonlat)

myfield = functionspace.create_field(name="swe_variables", variables=5, dtype = np.float64)
uvwh = atlas.make_view(myfield)

mycoords = functionspace.create_field(name="lonlat", variables=2, dtype = np.float64)
lonlat_sd = atlas.make_view(mycoords)

n_p = functionspace.size
search = Search(functionspace) #initializes the Search class with functionspace

z = xyz[:,2]
diff=np.zeros(n_p)
for i in range(n_p):

    value = np.cos(1/6)
    #z = xyz[:,2]
    diff[i] = np.abs(z[i] - value)

#print("Min:", np.argmin(diff), z[np.argmin(diff)])

index = np.argmin(diff)

nearest = search.nearest_indices_within_radius(index, myradius)

xyz_r = getneighcoords(nearest, xyz)

A, invA = constructA(xyz_r)


print("A:\n", A)
print("Determinant of A:", np.linalg.det(A))

Xj = xyz[index]
Dx, Dy, Dz = differentiation(invA,xyz_r,Xj)
fields = set_initial_conditions(xyz, n_p, ghost,lonlat)


#fill values in the functionspace variables
for n in range(n_p):
  if not ghost[n]:
    lonlat_sd[n,0] = lonlat[n,0]
    lonlat_sd[n,1] = lonlat[n,1]

    uvwh[n,0] = fields[n,0]
    uvwh[n,1] = fields[n,1]
    uvwh[n,2] = fields[n,2]
    uvwh[n,3] = fields[n,3]
    uvwh[n,4] = fields[n,4]

  else:
     lonlat_sd[n,0] = lonlat[n,1] = 0.0
     uvwh[n,0] = 0.0
     uvwh[n,1] = 0.0
     uvwh[n,2] = 0.0
     uvwh[n,3] = 0.0
     uvwh[n,4] = 0.0

myfield.halo_dirty = True
myfield.halo_exchange()

mycoords.halo_dirty = True
mycoords.halo_exchange()



## GRADIENT CHECKING####
#$\vec{t} = (-x,-y, \frac{x^2+y^2}{z})$

g = 1.5454e-6
h0 = 1000/6371000
R = 1.0/3
c_ana = np.zeros(3)
x = xyz[index,0]
y = xyz[index,1]
z = xyz[index,2]

uvwh_r = get_uvwh_r(uvwh, nearest)

#u = uvwh_r[:, 0]
#v = uvwh_r[:, 1]
#w = uvwh_r[:,2]
h = uvwh_r[:, 3]
           #rho = np.sqrt(x**2 +y**2 +((x**2+y**2)/z)**2)

print("h: ",h)

r = np.arccos(z)
vecnorth = np.array([-x, -y, ((x**2 +y**2)/z) ])
rho = np.linalg.norm(vecnorth)
vecnorth = vecnorth/rho

ana_gradient = g*(((h0*np.pi)/(2*R))*(-np.sin((np.pi*r)/R)))
                #print("yes")

c_ana[0] = vecnorth[0]*(ana_gradient)
c_ana[1] = vecnorth[1]*(ana_gradient)
c_ana[2] = vecnorth[2]*(ana_gradient)


termCx= (Dx)
termCy= (Dy)
termCz= (Dz)
    #ana_c =

term_c = np.row_stack([termCx,termCy,termCz])
termC = g*(np.matmul(term_c,h))

px, py, pz = getpxyz(Xj)


pxyz = np.row_stack([px,py,pz])

newtermc = -np.dot(pxyz,termC)
#termC[0] = -np.dot(px,termC[0])
#termC[1] = -np.dot(py,termC[1])
#termC[2] = -np.dot(pz,termC[2])

print("Numerical TermC:", newtermc)
print("Analytical TermC:", c_ana)
print("dot product:", np.dot(Xj,c_ana))
print("dot product numerical:", np.dot(Xj,newtermc))



atlas.finalize()
