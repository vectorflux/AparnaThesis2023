import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse
from all_functions import *

csr_mat=sparse.load_npz("icon_R2B4_wendland1_full_matrix.npz")
print("Non-zeros: ", csr_mat.nnz)
plt.spy(csr_mat,markersize=0.1)
plt.show()
invA=sparse.linalg.splu(csr_mat)
#sparse.save_npz("icon_R2B4_wendland1_full_inverse.npz", invA)
