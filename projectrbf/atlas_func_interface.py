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


def find_neighbors(myradius):#radius passed as argument

    atlas.initialize()

    #resolution = 100 * km / earth_radius
    resolution = myradius
    grid = atlas.Grid("H128")

    ###Create FunctionSpace
    functionspace = atlas.functionspace.PointCloud(grid, halo_radius=resolution * 2, geometry="UnitSphere")

    search = Search(functionspace)
    radius = resolution
    index = 0
    nearest = search.nearest_indices_within_radius(index, radius)
    if functionspace.part == 0:
        print("nearest global indices to local index", index, " ( global index", global_index[index], "): ",
              [global_index[n] for n in nearest])

    atlas.finalize()

    return nearest