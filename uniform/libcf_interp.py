import pycf
import numpy
import sys
from ctypes import byref, c_int, c_double, c_float, POINTER, c_char_p, c_void_p
import argparse
from functools import reduce

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

def createData(filename, prefix):

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
		                                     b"longitude", b"degrees_east", byref(lonCoordId))
	assert(ier == pycf.NC_NOERR)
	coordIds = (c_int * ndims)(latCoordId, lonCoordId)
	gridId = c_int()
	ier = pycf.nccf.nccf_def_grid(coordIds, (prefix + "grid").encode('UTF-8'), byref(gridId))
	assert(ier == pycf.NC_NOERR)

	dataId = c_int()
	read_data = 1
	ier = pycf.nccf.nccf_def_data_from_file(filename, gridId, b"air_temperature",
                                            read_data, byref(dataId))
	assert(ier == pycf.NC_NOERR)

	return gridId, dataId

def destroyData(dataId):
	gridId = c_int()
	ier = pycf.nccf.nccf_inq_data_gridid(dataId, byref(gridId))
	assert(ier == pycf.NC_NOERR)
	coordIds = (c_int * ndims)()
	ier = pycf.nccf.nccf_inq_grid_coordids(gridId, coordIds)
	assert(ier == pycf.NC_NOERR)

	ier = pycf.nccf.nccf_free_data(dataId)
	assert(ier == pycf.NC_NOERR)
	ier = pycf.nccf.nccf_free_grid(gridId)
	assert(ier == pycf.NC_NOERR)
	for i in range(ndims):
		ier = pycf.nccf.nccf_free_coord(coordIds[i])
		assert(ier == pycf.NC_NOERR)


srcGridId, srcDataId = createData(src_file, b"src")
dstGridId, dstDataId = createData(dst_file, b"dst")

# compute the interpolation weights
regridId = c_int()
ier = pycf.nccf.nccf_def_regrid(srcGridId, dstGridId, byref(regridId))
assert(ier == pycf.NC_NOERR)
nitermax = 100
tolpos = c_double(1.e-10)
ier = pycf.nccf.nccf_compute_regrid_weights(regridId,
                                            nitermax, tolpos)
assert(ier == pycf.NC_NOERR)
ntargets = c_int()
ier = pycf.nccf.nccf_inq_regrid_ntargets(regridId, byref(ntargets))
assert(ier == pycf.NC_NOERR)
nvalid = c_int()
ier = pycf.nccf.nccf_inq_regrid_nvalid(regridId, byref(nvalid))
assert(ier == pycf.NC_NOERR)

# store the reference data values
xtypep = c_int()
dstDataPtr = POINTER(c_double)()
fillValuePtr = c_void_p()
ier = pycf.nccf.nccf_get_data_pointer(dstDataId, byref(xtypep),
                                      byref(dstDataPtr), byref(fillValuePtr))
assert(ier == pycf.NC_NOERR)
dims = (c_int * ndims)()
ier = pycf.nccf.nccf_inq_data_dims(dstDataId, dims)
assert(ier == pycf.NC_NOERR)
ntot = reduce(lambda x, y: x*y, dims[:], 1)
dstDataRef = numpy.zeros((ntot,), numpy.float64)
for i in range(ntot): 
	dstDataRef[i] = dstDataPtr[i] # copy

# interpolate (CRASHES!)
ier = pycf.nccf.nccf_apply_regrid(regridId, srcDataId, dstDataId)
assert(ier == pycf.NC_NOERR)

# compute error
error = 0.0
for i in range(ntot):
	error += abs(dstDataPtr[i] - dstDataRef[i])
error /= float(ntot)
print('interpolation error: {} (# targets = {}, # valid points = {})'.format(error, 
	  ntargets.value, nvalid.value))

# clean up
destroyData(srcDataId)
destroyData(dstDataId)
