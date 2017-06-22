import iris
import ants
import ants.tests.stock as stock
from matplotlib import pylab

dlon, dlat = 10., 20.

cube = stock.geodetic(shape=(10, 20), with_bounds=True, circular=True, xlim=(-180., 180.), ylim=(-90., 90.))
#cube.coord_system = iris.coord_systems.RotatedGeogCS(grid_north_pole_latitude=20.0, grid_north_pole_longitude=30.0)


lonMid2d, latMid2d = iris.analysis.cartography.get_xy_grids(cube)
lonCrn2d, latCrn2d = iris.analysis.cartography.get_xy_contiguous_bounded_grids(cube)

lonCrnRotated, latCrnRotated = iris.analysis.cartography.rotate_pole(lonCrn2d, latCrn2d,
	                                                                 pole_lon=180.0 + dlon, pole_lat=90.0 + dlat)

m, n = lonCrn2d.shape
pylab.figure(1)
for j in range(m):
	pylab.plot(lonCrnRotated[j, :], latCrnRotated[j, :], 'g-')
for i in range(n):
	pylab.plot(lonCrnRotated[:, i], latCrnRotated[:, i], 'r-')
pylab.savefig('latlonRotated.png')

pylab.figure(2)
for j in range(m):
	pylab.plot(lonCrn2d[j, :], latCrn2d[j, :], 'g-')
for i in range(n):
	pylab.plot(lonCrn2d[:, i], latCrn2d[:, i], 'r-')
pylab.savefig('latlon.png')

print(cube.coord_system())
#iris.save([cube], 'cube.nc')
