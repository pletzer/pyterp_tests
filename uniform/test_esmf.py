"""
example how to create a structured grid from lat-lon arrays
"""

import numpy
import ESMF
lats = numpy.linspace(-90., 90., 11)
lons = numpy.linspace(0., 360., 21)

cellDims = numpy.array([len(lats) - 1, len(lons) - 1])
grid = ESMF.Grid(max_index=cellDims)

latIndex, lonIndex = 0, 1
grid.add_coords(staggerloc=ESMF.StaggerLoc.CORNER, coord_dim=latIndex)
grid.add_coords(staggerloc=ESMF.StaggerLoc.CORNER, coord_dim=lonIndex)

coordLat = grid.get_coords(coord_dim=latIndex, staggerloc=ESMF.StaggerLoc.CORNER)
coordLon = grid.get_coords(coord_dim=lonIndex, staggerloc=ESMF.StaggerLoc.CORNER)

iBegLat = grid.lower_bounds[ESMF.StaggerLoc.CORNER][latIndex]
iEndLat = grid.upper_bounds[ESMF.StaggerLoc.CORNER][latIndex]
iBegLon = grid.lower_bounds[ESMF.StaggerLoc.CORNER][lonIndex]
iEndLon = grid.upper_bounds[ESMF.StaggerLoc.CORNER][lonIndex]

for j in range(iBegLat, iEndLat):
	for i in range(iBegLon, iEndLon):
		coordLat[j, i] = lats[j]
		coordLon[j, i] = lons[i]

print(coordLat)
print(coordLon)

field = ESMF.Field(grid, name='air_temperature', staggerloc=ESMF.StaggerLoc.CORNER)
field.data[...] = numpy.sin(numpy.pi * coordLat/180.)

print(field.data.shape)