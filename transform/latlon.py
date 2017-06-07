import iris
import ants
import ants.tests.stock as stock
from matplotlib import pylab

cube = stock.geodetic(shape=(10, 20), with_bounds=True, circular=True, xlim=(-180., 180.), ylim=(-90., 90.))
#cube.coord_system = iris.coord_systems.RotatedGeogCS(grid_north_pole_latitude=20.0, grid_north_pole_longitude=30.0)


lonMid2d, latMid2d = iris.analysis.cartography.get_xy_grids(cube)
lonCrn2d, latCrn2d = iris.analysis.cartography.get_xy_contiguous_bounded_grids(cube)

lonCrnRotated, latCrnRotated = iris.analysis.cartography.rotate_pole(lonCrn2d, latCrn2d, pole_lon=30.0, pole_lat=20.0)

m, n = lonCrn2d.shape
for j in range(m):
	pylab.plot(lonCrnRotated[j, :], latCrnRotated[j, :], 'g-')
for i in range(n):
	pylab.plot(lonCrnRotated[:, i], latCrnRotated[:, i], 'r-')
pylab.savefig('latlon.png')

print(cube.coord_system())
#iris.save([cube], 'cube.nc')
