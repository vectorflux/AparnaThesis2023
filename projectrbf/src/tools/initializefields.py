import atlas4py as atlas

class fields:
    def __init__(self, functionspace):
        self.functionspace = functionspace
        self.lonlat = atlas.make_view(self.functionspace.lonlat)
        self.fields = atlas.make_view(self.functionspace.fields)

        #u, v, w, h



def init_points(Dx,Dy,Dz):




