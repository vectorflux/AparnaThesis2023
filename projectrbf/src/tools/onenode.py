#### Just for one node at 1717 #####


#!/usr/bin/env python3
import sys
import atlas4py as atlas
import numpy as np
import time


from scipy import linalg

import pandas as pd

from operator_matrices import *
from evaluate_wendland1 import *
from atlas_func_interface import *
from read_netcdf_file import *
from A_matrices import *
from initializefields import *
from visualization import *
#from rk4 import *
from construct_rhsd import *

np.set_printoptions(threshold=sys.maxsize)

atlas.initialize() #initializes atlas and MPI


myradius = 0.100

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

myfield = functionspace.create_field(name="swe_variables", variables=4, dtype = np.float64)
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

index = 1717

print("Node coordinates:\n", xyz[index,:])

nearest = search.nearest_indices_within_radius(index, myradius)

xyz_r = getneighcoords(nearest, xyz)


print("Nearest neighbors:\n",xyz_r)

A, invA = constructA(xyz_r,myradius)


print("Maximum of A:", np.max(A))
#print("Len of A:", len(A))
#for i in range(len(A)):
    #for j in range(len(A)):
        #print(A[i][j])
#print("A:\n", A)
#print("Determinant of A:", np.linalg.det(A))


DF = pd.DataFrame(A)

DF.to_csv("Amatrix.csv")


print("A:\n", A)
print("Determinant of A:", np.linalg.det(A))
print("max(A): ", A.max() )


Xj = xyz[index]
Dx, Dy, Dz = differentiation(invA,xyz_r,Xj,myradius)
fields, f = set_initial_conditions(xyz, n_p, ghost,lonlat)


#fill values in the functionspace variables
for n in range(n_p):
  if not ghost[n]:
    lonlat_sd[n,0] = lonlat[n,0]
    lonlat_sd[n,1] = lonlat[n,1]

    uvwh[n,0] = fields[n,0]
    uvwh[n,1] = fields[n,1]
    uvwh[n,2] = fields[n,2]
    uvwh[n,3] = fields[n,3]

  else:
     lonlat_sd[n,0] = lonlat[n,1] = 0.0
     uvwh[n,0] = 0.0
     uvwh[n,1] = 0.0
     uvwh[n,2] = 0.0
     uvwh[n,3] = 0.0

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

u = uvwh_r[:, 0]
#v = uvwh_r[:, 1]
#w = uvwh_r[:,2]
h = uvwh_r[:, 3]
           #rho = np.sqrt(x**2 +y**2 +((x**2+y**2)/z)**2)


print("h: ",h)
c = np.matmul(invA,h)
print("c = invA*h: ",c)
h_approx = evaluate_wendland1(c,xyz_r,Xj,myradius)
print("h approximation: ", h_approx )
print("h true value: ", uvwh[index,3] )
print("h relerr:", abs(h_approx - uvwh[index,3])/abs(uvwh[index,3]) )

gradient_j = gradient(c,xyz_r,Xj,myradius)
print("Unprojected gradient approx:\n", gradient(c,xyz_r,Xj,myradius) )

proj_gradient_j = gradient_j - Xj*np.dot(gradient_j,Xj)
print("Projected gradient approx 1:\n", g*proj_gradient_j )
print("Should be zero: ", np.dot(proj_gradient_j,Xj) )

# Not sure why this is in the order of 1.e-5 --> too large
#print("Node j: ", Xj )
#print("Component in zonal direction (should be small): ", np.dot( gradient_j, [ -Xj[2], Xj[1], 0 ] ) )

#proj_gradient(c,xyz_r,Xj,myradius)
#proj_gradient_j = proj_gradient(c,xyz_r,Xj,myradius)
#print("Projected gradient approx 2:\n", proj_gradient_j )
#print("Is not zero: ", np.dot(proj_gradient_j,Xj) )

#print("h: ",h)


r = np.arccos(z)
vecnorth = np.array([-x, -y, ((x**2 +y**2)/z) ])
rho = np.linalg.norm(vecnorth)
vecnorth = vecnorth/rho


ana_gradient = ((h0*np.pi)/(2*R))*(-np.sin((np.pi*r)/R))


c_ana[0] = vecnorth[0]*(ana_gradient)
c_ana[1] = vecnorth[1]*(ana_gradient)
c_ana[2] = vecnorth[2]*(ana_gradient)


c_ana    = g*c_ana   # create analytic TermC


termCx= (Dx)
termCy= (Dy)
termCz= (Dz)
    #ana_c =

term_c = np.row_stack([termCx,termCy,termCz])

termC = g*(np.matmul(term_c,h))
#termC = g*gradient_j
px, py, pz = getpxyz(Xj)
#print("pxyz:", px,py,pz)

pxyz = np.row_stack([px,py,pz])
newtermc = np.dot(pxyz,termC)/myradius   # normalize here, since normalization factor not in matrix
#newtermc = termC - Xj*np.dot(termC,Xj)

##termC[0] = -np.dot(px,termC[0])
#termC[1] = -np.dot(py,termC[1])
#termC[2] = -np.dot(pz,termC[2])

print("Numerical TermC:", newtermc)
print("Analytical TermC:", c_ana)
relerr_gradient = np.linalg.norm(newtermc/g - c_ana/g)/np.linalg.norm(c_ana/g)
print("relative difference in gradient", relerr_gradient)

print("arc deviation of gradient", np.arccos(np.dot(newtermc/g,c_ana/g) / (np.linalg.norm(newtermc/g) * np.linalg.norm(c_ana/g)) ) )

print("dot product:", np.dot(Xj,c_ana))
print("dot product numerical:", np.dot(Xj,newtermc))

newtermc1 = np.dot(pxyz,termC)
newtermc2 = termC - Xj*np.dot(termC,Xj)

#termC[0] = -np.dot(px,termC[0])
#termC[1] = -np.dot(py,termC[1])
#termC[2] = -np.dot(pz,termC[2])

#print("Numerical TermC with pxyz:", newtermc1/myradius)
#print("Numerical TermC correct:", newtermc2/myradius)
#print("Analytical TermC:", c_ana)
##print("dot product:", np.dot(Xj,c_ana))
#print("dot product of projectors:", np.dot(Xj,pxyz))


#print("dot product numerical:", np.dot(Xj,newtermc))


print("neighborhood: ", np.shape(xyz_r))

atlas.finalize()
