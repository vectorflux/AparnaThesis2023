import numpy as np
from scipy import sparse
from all_functions import *

myradius = 0.100 # radius of neighborhood on unit sphere
#r = 0.850 # radius of neighborhood on unit sphere
coords = read_data_netcdf() #Read data from a netcdf file

coord_size = len(coords)

data = []
indices = []
indptr = [0]

#Main Loop which iterates over all points and calls different functions
for id in range(coord_size):
    my_list = create_neighborhood(myradius, coords, id) #Myimplementation

    #print('my list', my_list)
    print("row:", id, "neighbors:", len(my_list))
    #print(coords[my_list[id]])

    add_row( data, indices, indptr, my_list, getrow(id, my_list, coords, myradius) )
    print('***End of loop ', id, 'nonzeros ', len(my_list), '***')

csr_mat = sparse.csr_matrix((data, indices, indptr), shape=(coord_size, coord_size))
print('*** CSR matrix now constructed ***')

csr_mat = csr_mat.asfptype()
print('csr_mat', csr_mat)
sparse.save_npz("icon_R2B4_wendland1_full_matrix.npz", csr_mat)

ew1, ev = sparse.linalg.eigs(csr_mat, which='LM')
print('max(eigs)', ew1)
print('*** Calculated largest eigenvalues by magnitude ***')
ew2, ev = sparse.linalg.eigs(csr_mat, sigma=1e-8)
print('*** Calculated smallest eigenvalues by magnitude ***')
cond_A = max(abs(ew1)) / min(abs(ew2))
print('min(eigs)', ew2)
print('***Condition of A ', cond_A, '***')

invA=sparse.linalg.splu(csr_mat)
sparse.save_npz("icon_R2B4_wendland1_full_inverse.npz", invA)


