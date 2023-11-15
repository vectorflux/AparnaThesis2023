import atlas4py as atlas
import math
import numpy as np

#atlas.initialize() #initializes atlas and MPI

# constants
earth_radius = 6371229.
km = 1000.

class Search:
    def __init__(self, functionspace):
        self.functionspace = functionspace
        self.lonlat = atlas.make_view(self.functionspace.lonlat)
        self.kdtree = atlas.IndexKDTree(geometry="UnitSphere")
        self.kdtree.build(self.lonlat)

    def nearest_indices_within_radius(self, i, radius):
        # radius is in metres on Earth geometry
        closest_indices = self.kdtree.closest_indices_within_radius(lon=self.lonlat[i, 0], lat=self.lonlat[i, 1],
                                                                    radius=radius)
        return closest_indices[1:]  # first point is "i" itself, so skip it



#def update_fields(self,):

