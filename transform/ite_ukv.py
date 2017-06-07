import iris
import ants
import ants.tests.stock as stock
from matplotlib import pylab

cube = iris.load_cube('ite_ukv.nc')

# get the modified lats and lons
grid_lat = cube.coords('grid_latitude')
grid_lon = cube.coords('grid_longitude')

# get the coordinate reference system
crs = cube.coord_system()

# identify the type of coordinate system as being rotated pole
dlon, dlat = 0., 0.
if isinstance(crs, iris.coord_systems.RotatedGeogCS):
	dlat = crs.grid_north_pole_latitude
	dlon = crs.grid_north_pole_longitude

# get the rotated pole coordinates
grid_lon2d, grid_lat2d = iris.analysis.cartography.get_xy_contiguous_bounded_grids(cube)
# compute the true lat/lon
lon2d, lat2d = iris.analysis.cartography.rotate_pole(grid_lon2d, grid_lat2d, 
	                                                 pole_lon=dlon, pole_lat=dlat)

m, n = lon2d.shape
pylab.figure(0)
for j in range(m):
	pylab.plot(lon2d[j, :], lat2d[j, :], 'g-')
for i in range(n):
	pylab.plot(lon2d[:, i], lat2d[:, i], 'r-')
pylab.savefig('true_latlon.png')

pylab.figure(1)
for j in range(m):
	pylab.plot(grid_lon2d[j, :], grid_lat2d[j, :], 'g-')
for i in range(n):
	pylab.plot(grid_lon2d[:, i], grid_lat2d[:, i], 'r-')
pylab.savefig('grid_latlon.png')


print(grid_lat)
print(grid_lon)
