from __future__ import print_function
import ESMF
import iris
import numpy
import sys
from ctypes import byref, c_int, c_double, c_float, POINTER, c_char_p, c_void_p
import argparse
from functools import reduce
import time

parser = argparse.ArgumentParser(description='Interpolate using libcf')
parser.add_argument('--src_file', type=str, dest='src_file', default='src.nc',
                    help='Source data file name')
parser.add_argument('--dst_file', type=str, dest='dst_file', default='dst.nc',
                    help='Destination data file name')

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
	cube = iris.load_cube(filename, 'air_temperature')
	coords = cube.coords()
	lats = coords[0].points
	lons = coords[1].points
	
	# NOTE fortran ordering here

	# create the ESMF grid object

	latIndex, lonIndex = 0, 1 # or 1, 0????
	cellDims = numpy.array([len(lats) - 1, len(lons) - 1])
	grid = ESMF.Grid(max_index=cellDims)

	grid.add_coords(staggerloc=ESMF.StaggerLoc.CORNER, coord_dim=latIndex)
	grid.add_coords(staggerloc=ESMF.StaggerLoc.CORNER, coord_dim=lonIndex)

	coordLat = grid.get_coords(coord_dim=latIndex, staggerloc=ESMF.StaggerLoc.CORNER)
	coordLon = grid.get_coords(coord_dim=lonIndex, staggerloc=ESMF.StaggerLoc.CORNER)

	# get the local start/end index sets
	iBegLat = grid.lower_bounds[ESMF.StaggerLoc.CORNER][latIndex]
	iEndLat = grid.upper_bounds[ESMF.StaggerLoc.CORNER][latIndex]
	iBegLon = grid.lower_bounds[ESMF.StaggerLoc.CORNER][lonIndex]
	iEndLon = grid.upper_bounds[ESMF.StaggerLoc.CORNER][lonIndex]

	# set the coordinates
	for j in range(iBegLat, iEndLat):
		for i in range(iBegLon, iEndLon):
			coordLat[j, i] = lats[j]
			coordLon[j, i] = lons[i]

	# create field
	field = ESMF.Field(grid, name="air_temperature", 
		               staggerloc=ESMF.StaggerLoc.CORNER)
	print(field.data.shape)
	field.data[...] = cube.data[:]

	nodeDims = (iEndLat - iBegLat, iEndLon - iBegLon)

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
	                            regrid_method=None, pole_method=None, regrid_pole_npoints=None, 
	                            line_type=None, norm_type=None, unmapped_action=None, 
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
print('emsf interpolation:')
print('\tsrc: {} ntot: {}'.format(srcNodeDims, srcNtot))
print('\tdst: {} ntot: {}'.format(dstNodeDims, dstNtot))
print('interpolation error: {}'.format(error))
print('time stats:')
for k, v in timeStats.items():
	print('\t{0:<32} {1:>.3g} sec'.format(k, v))
print()

# clean up
# nothing to do