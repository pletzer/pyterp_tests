from __future__ import print_function
import ESMF
import iris
import numpy
import sys
from ctypes import byref, c_int, c_double, c_float, POINTER, c_char_p, c_void_p
import argparse
from functools import reduce
import time

LAT_INDEX, LON_INDEX = 1, 0

parser = argparse.ArgumentParser(description='Interpolate using ESMF')
parser.add_argument('--src_file', type=str, dest='src_file', default='src.nc',
                    help='Source data file name')
parser.add_argument('--dst_file', type=str, dest='dst_file', default='dst.nc',
                    help='Destination data file name')
parser.add_argument('--plot', dest='plot', action='store_true', help='Plot')

args = parser.parse_args()

if args.src_file is '':
    print('ERROR: must provide source data file name')
    parser.print_help()
    sys.exit(1)

if args.dst_file is '':
    print('ERROR: must provide destination data file name')
    parser.print_help()
    sys.exit(1)

src_file = args.src_file.encode('UTF-8') # python3
dst_file = args.dst_file.encode('UTF-8') # python3
ndims = 2

def createData(filename, prefix):
    # use iris to read in the data
    # then pass the array to the ESMF API
    pointCube = iris.load(filename, iris.Constraint(cube_func = lambda c: c.var_name == 'pointData'))[0]
    cellCube = iris.load(filename, iris.Constraint(cube_func = lambda c: c.var_name == 'cellData'))[0]
    coords = pointCube.coords()
    lats = coords[0].points
    lons = coords[1].points
    
    # NOTE fortran ordering here

    # create the ESMF grid object

    cellDims = numpy.array([len(lons) - 1, len(lats) - 1])
    grid = ESMF.Grid(max_index=cellDims)

    grid.add_coords(staggerloc=ESMF.StaggerLoc.CORNER, coord_dim=LAT_INDEX)
    grid.add_coords(staggerloc=ESMF.StaggerLoc.CORNER, coord_dim=LON_INDEX)

    coordLat = grid.get_coords(coord_dim=LAT_INDEX, staggerloc=ESMF.StaggerLoc.CORNER)
    coordLon = grid.get_coords(coord_dim=LON_INDEX, staggerloc=ESMF.StaggerLoc.CORNER)

    # get the local start/end index sets
    iBegLat = grid.lower_bounds[ESMF.StaggerLoc.CORNER][LAT_INDEX]
    iEndLat = grid.upper_bounds[ESMF.StaggerLoc.CORNER][LAT_INDEX]
    iBegLon = grid.lower_bounds[ESMF.StaggerLoc.CORNER][LON_INDEX]
    iEndLon = grid.upper_bounds[ESMF.StaggerLoc.CORNER][LON_INDEX]

    # set the coordinates
    coordLat[...] = numpy.outer(numpy.ones((iEndLon - iBegLon,), coordLon.dtype), lats[iBegLat:iEndLat])
    coordLon[...] = numpy.outer(lons[iBegLon:iEndLon], numpy.ones((iEndLat - iBegLat,), coordLon.dtype))

    # create field
    field = ESMF.Field(grid, name="air_temperature", 
                   staggerloc=ESMF.StaggerLoc.CENTER)
    field.data[...] = cellCube.data.transpose()

    nodeDims = (iEndLon - iBegLon, iEndLat - iBegLat)

    return grid, field, nodeDims

timeStats = {
    'weights': float('nan'),
    'evaluation': float('nan'),
}

srcGrid, srcData, srcNodeDims = createData(src_file, b"src")
dstGrid, dstData, dstNodeDims = createData(dst_file, b"dst")

# save the reference (exact) field data
dstDataRef = dstData.data.copy()

# compute the interpolation weights
tic = time.time()
regrid = ESMF.api.regrid.Regrid(srcData, dstData,
                                src_mask_values=None, dst_mask_values=None,
                                regrid_method=ESMF.api.constants.RegridMethod.CONSERVE, 
                                pole_method=None, regrid_pole_npoints=None, 
                                line_type=None, norm_type=None, 
                                unmapped_action=ESMF.api.constants.UnmappedAction.IGNORE, 
                                ignore_degenerate=None, src_frac_field=None, dst_frac_field=None)
timeStats['weights'] = time.time() - tic

# interpolate
tic = time.time()
regrid(srcData, dstData)
timeStats['evaluation'] = time.time() - tic

# compute error
srcNtot = len(srcData.data.flat)
dstNtot = len(dstData.data.flat)
error =  numpy.sum(abs(dstData.data - dstDataRef)) / float(dstNtot)
print('emsf conservative interpolation:')
print('\tsrc: {} ntot: {}'.format(srcNodeDims, srcNtot))
print('\tdst: {} ntot: {}'.format(dstNodeDims, dstNtot))
print('interpolation error: {:.3g}'.format(error))
totTime = 0.0
print('time stats:')
for k, v in timeStats.items():
    print('\t{0:<32} {1:>.3g} sec'.format(k, v))
    totTime += v
print('\t{0:<32} {1:>.3g} sec'.format('total', totTime))

# check sum
checksum = numpy.sum(dstData.data, axis=None)
print('check sum: {:.15g}'.format(checksum))

# plot
if args.plot:
    lats = dstGrid.get_coords(coord_dim=LAT_INDEX, staggerloc=ESMF.StaggerLoc.CORNER)
    lons = dstGrid.get_coords(coord_dim=LON_INDEX, staggerloc=ESMF.StaggerLoc.CORNER)
    from matplotlib import pylab
    pylab.pcolor(lons, lats, dstData.data)
    pylab.show()

# clean up
# nothing to do
