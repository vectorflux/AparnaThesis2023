import atlas4py as atlas
import math
import numpy as np

atlas.initialize() #initializes atlas and MPI

# constants
earth_radius = 6371229.
km = 1000.

class Search:
    def __init__(self, functionspace):
        self.functionspace = functionspace
        self.lonlat = atlas.make_view(self.functionspace.lonlat)
        self.kdtree = atlas.IndexKDTree(geometry="UnitSphere")
        self.kdtree.build(lonlat)

    def nearest_indices_within_radius(self, i, radius):
        # radius is in metres on Earth geometry
        closest_indices = self.kdtree.closest_indices_within_radius(lon=self.lonlat[i, 0], lat=self.lonlat[i, 1],
                                                                    radius=radius)
        return closest_indices[1:]  # first point is "i" itself, so skip it



class fields:
    def __init__(self, functionspace):
        self.functionspace = functionspace
        self.lonlat = atlas.make_view(self.functionspace.lonlat)
        self.fields = atlas.make_view(self.functionspace.fields)

        #u, v, w, h




def find_neighbors(myradius,functionspace):#radius passed as argument

    n_p = functionspace.size
    search = Search(functionspace)

    for id in range(functionspace.size): #n_p
        if ghost[id] == 0:
            nearest = search.nearest_indices_within_radius(id, myradius)
            #print("nearest global indices to local index", id, " ( global index", global_index[id], "): ",
            #  [global_index[n] for n in nearest])


    return nearest

#def update_fields(self,):

