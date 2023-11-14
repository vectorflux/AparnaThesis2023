import atlas4py as atlas
from operator_matrices import *

class fields:
    def __init__(self, functionspace):
        self.functionspace = functionspace
        self.lonlat = atlas.make_view(self.functionspace.lonlat)
        self.fields = atlas.make_view(self.functionspace.fields)

        #u, v, w, h


def initialize_fields(uvwh):




def construct_rhs(uvwh, allD,xyz,nrj_size_list):

    px, py, pz = projection(xyz_r)

    Ru =

    Rv =

    Rz =


