import pycf
import numpy
import sys
from ctypes import byref, c_int, c_double, c_float, POINTER, c_char_p
import argparse

parser = argparse.ArgumentParser(description='Interpolate using libcf')
parser.add_argument('--src_file', type=str, dest='src_file', default='',
                    help='Source data file name')
parser.add_argument('--dst_file', type=str, dest='dst_file', default='',
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

def createUniformData(filename, prefix):

	latAxisId = c_int()
	lonAxisId = c_int()
	latCoordId = c_int()
	lonCoordId = c_int()
	gridId = c_int()
	dataId = c_int()

	ier = pycf.nccf.nccf_def_axis_from_file(filename, b"latitude", byref(latAxisId))
	assert(ier == pycf.NC_NOERR)
	ier = pycf.nccf.nccf_def_axis_from_file(filename, b"longitude", byref(lonAxisId))
	assert(ier == pycf.NC_NOERR)

	axisIds = (c_int * ndims)(latAxisId, lonAxisId)
	ier = pycf.nccf.nccf_def_coord_from_axes(ndims, axisIds, 0, (prefix + "_lats").encode('UTF-8'),
		                                     b"latitude", b"degrees_north", byref(latCoordId))
	assert(ier == pycf.NC_NOERR)
	ier = pycf.nccf.nccf_def_coord_from_axes(ndims, axisIds, 1, (prefix + "_lons").encode('UTF-8'), 
		                                     "longitude", "degrees_east", byref(lonCoordId))
	assert(ier == pycf.NC_NOERR)
	coordIds = (c_int * ndims)(latCoordId, lonCoordId)
	gridId = c_int()
	ier = pycf.nccf.nccf_def_grid(coordIds, "srcGrid", byref(gridId))
	assert(ier == pycf.NC_NOERR)

	return gridId

srcGridId = createUniformData(src_file, "src")
dstGridId = createUniformData(src_file, "dst")

# interpolate
regridId = c_int()
ier = pycf.nccf.nccf_def_regrid(srcGridId, dstGridId, byref(regridId))
assert(ier == pycf.NC_NOERR)
nitermax = 100
tolpos = c_double(1.e-10)
ier = pycf.nccf.nccf_compute_regrid_weights(regridId,
                                            nitermax, tolpos)
assert(ier == pycf.NC_NOERR)

# compute error

# clean up