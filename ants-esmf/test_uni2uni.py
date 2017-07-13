import sys
import iris
import ants
from ants.regrid import ConservativeESMF
from matplotlib import pylab
import numpy

def get_cell_areas(cube):
	lon_bnds = cube.coord('longitude').bounds
	lat_bnds = cube.coord('latitude').bounds
	dlon = lon_bnds[:, 1] - lon_bnds[:, 0]
	dlat = lat_bnds[:, 1] - lat_bnds[:, 0]
	return numpy.outer(dlat, dlon)

src_cube = iris.load('uni.nc', iris.Constraint(cube_func = lambda c: c.var_name == 'cellData'))[0]
tgt_cube = iris.load('uni2.nc', iris.Constraint(cube_func = lambda c: c.var_name == 'cellData'))[0]
rgrd_factory = ConservativeESMF()
rgrd = rgrd_factory.regridder(src_cube, tgt_cube)

# take square of sine to make it easier to take sums
src_cube.data = src_cube.data**2
tgt_cube.data = tgt_cube.data**2

# regrid
res_cube = rgrd(src_cube)
print res_cube
iris.save(res_cube, 'uni_uni2.nc')

# run some checks
src_sum = (src_cube.data * get_cell_areas(src_cube)).sum()
res_sum = (res_cube.data * get_cell_areas(res_cube)).sum()
print 'src_sum = {} res_sum = {}'.format(src_sum, res_sum)

# plot
for s, cube in [('src', src_cube), ('tgt', tgt_cube), ('res', res_cube)]:
	print s
	pylab.figure()
	lon = cube.coord('longitude').points
	lat = cube.coord('latitude').points
	pylab.pcolor(lon, lat, cube.data)
	pylab.title(s)

pylab.show()

# some diffs are expected because of the curvature of the grids
#assert(abs(src_sum - res_sum) < 1.e-6)