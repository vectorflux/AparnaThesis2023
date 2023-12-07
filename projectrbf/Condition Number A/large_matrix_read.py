import numpy as np
from scipy import sparse
from all_functions import *

csr_mat=sparse.load_npz("icon_R2B4_wendland1_full_matrix.npz")

invA=sparse.linalg.splu(csr_mat)
sparse.save_npz("icon_R2B4_wendland1_full_inverse.npz", invA)
