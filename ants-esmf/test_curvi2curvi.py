import sys
import iris
import ants
from ants.regrid import ConservativeESMF
from matplotlib import pylab
import numpy

def get_cell_areas(cube):
    """
    areas of the spehrical quad based on the two spherical triangles
    http://mathworld.wolfram.com/SphericalTriangle.html
    """
    lon_bnd = cube.coord('longitude').bounds
    lat_bnd = cube.coord('latitude').bounds
    if len(lon_bnd.shape) == 2:
        # 1d axes, make it 2d
        nlat, nlon = lat_bnd.shape[0], lon_bnd.shape[0]

        lat2d_bnd = numpy.zeros((nlat, nlon, 4))
        lat2d_bnd[:, :, 0] = numpy.outer(lat_bnd[:, 0], numpy.ones((nlon,), numpy.float64))
        lat2d_bnd[:, :, 1] = numpy.outer(lat_bnd[:, 0], numpy.ones((nlon,), numpy.float64))
        lat2d_bnd[:, :, 2] = numpy.outer(lat_bnd[:, 1], numpy.ones((nlon,), numpy.float64))
        lat2d_bnd[:, :, 3] = numpy.outer(lat_bnd[:, 1], numpy.ones((nlon,), numpy.float64))

        lon2d_bnd = numpy.zeros((nlat, nlon, 4))
        lon2d_bnd[:, :, 0] = numpy.outer(numpy.ones((nlat,), numpy.float64), lon_bnd[:, 0])
        lon2d_bnd[:, :, 1] = numpy.outer(numpy.ones((nlat,), numpy.float64), lon_bnd[:, 1])
        lon2d_bnd[:, :, 2] = numpy.outer(numpy.ones((nlat,), numpy.float64), lon_bnd[:, 1])
        lon2d_bnd[:, :, 3] = numpy.outer(numpy.ones((nlat,), numpy.float64), lon_bnd[:, 0])

    else:
        # already 2d
        lat2d_bnd = lat_bnd
        lon2d_bnd = lon_bnd
        nlat, nlon = lat2d_bnd.shape[:2]

    # cartesian positions
    pos_bnd = numpy.zeros((nlat, nlon, 4, 3), numpy.float64)
    # x
    pos_bnd[..., 0] = numpy.cos(lat2d_bnd*numpy.pi/180.) * numpy.cos(lon2d_bnd*numpy.pi/180.)
    # y
    pos_bnd[..., 1] = numpy.cos(lat2d_bnd*numpy.pi/180.) * numpy.sin(lon2d_bnd*numpy.pi/180.)
    # z
    pos_bnd[..., 2] = numpy.sin(lat2d_bnd*numpy.pi/180.)

    # compute the angles for each corner
    angles_bnd = numpy.zeros((nlat, nlon, 4), numpy.float64)
    angles_bnd[:, :, 0] = numpy.arccos(((pos_bnd[..., 1, :] - pos_bnd[..., 0, :])*(pos_bnd[..., 3, :] - pos_bnd[..., 0, :])).sum(axis=2))
    angles_bnd[:, :, 1] = numpy.arccos(((pos_bnd[..., 2, :] - pos_bnd[..., 1, :])*(pos_bnd[..., 0, :] - pos_bnd[..., 1, :])).sum(axis=2))
    angles_bnd[:, :, 2] = numpy.arccos(((pos_bnd[..., 3, :] - pos_bnd[..., 2, :])*(pos_bnd[..., 1, :] - pos_bnd[..., 2, :])).sum(axis=2))
    angles_bnd[:, :, 3] = numpy.arccos(((pos_bnd[..., 0, :] - pos_bnd[..., 3, :])*(pos_bnd[..., 2, :] - pos_bnd[..., 3, :])).sum(axis=2))

    # for the quad it's adding all the angles - 2*pi
    spherical_excess = angles_bnd.sum(axis=2) - 2*numpy.pi
    # assuming radius = 1
    return spherical_excess


src_cube = iris.load('curvi.nc', iris.Constraint(cube_func = lambda c: c.var_name == 'cellData'))[0]
tgt_cube = iris.load('curvi2.nc', iris.Constraint(cube_func = lambda c: c.var_name == 'cellData'))[0]
rgrd_factory = ConservativeESMF()
rgrd = rgrd_factory.regridder(src_cube, tgt_cube)

# take square of sine to make it easier to take sums
src_cube.data = src_cube.data**2
tgt_cube.data = tgt_cube.data**2

# regrid
res_cube = rgrd(src_cube)
print res_cube
iris.save(res_cube, 'curvi_curvi.nc')

# run some checks
src_areas = get_cell_areas(src_cube)
res_areas = get_cell_areas(res_cube)
print 'src_area = {} res_area = {}'.format(src_areas.sum(), res_areas.sum())
src_sum = (src_cube.data * src_areas).sum()
res_sum = (res_cube.data * res_areas).sum()
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