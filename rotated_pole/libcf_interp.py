from __future__ import print_function
import pycf
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
parser.add_argument('--tolpos', type=float, dest='tolpos', default=1.e-8,
	                help='Tolerance in target space')
parser.add_argument('--nitermax', type=int, dest='nitermax', default=100,
	                help='Max number of iterations')

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

	latCoordId = c_int()
	lonCoordId = c_int()
	gridId = c_int()
	dataId = c_int()

	ier = pycf.nccf.nccf_def_coord_from_file(filename, 
		                                     b"latitude",
		                                     byref(latCoordId))
	assert(ier == pycf.NC_NOERR)
	ier = pycf.nccf.nccf_def_coord_from_file(filename, 
		                                     b"longitude",
		                                     byref(lonCoordId))
	assert(ier == pycf.NC_NOERR)

	coordIds = (c_int * ndims)(latCoordId, lonCoordId)
	gridId = c_int()
	ier = pycf.nccf.nccf_def_grid(coordIds, prefix + b"grid", byref(gridId))
	assert(ier == pycf.NC_NOERR)

	periodicity_lengths = (c_double * ndims)()
	ier = pycf.nccf.nccf_inq_grid_periodicity(gridId, periodicity_lengths)
	assert(ier == pycf.NC_NOERR)
	print('periodicity lengths: {}'.format(periodicity_lengths[:]))

	dataId = c_int()
	read_data = 1
	ier = pycf.nccf.nccf_def_data_from_file(filename, gridId, b"air_temperature",
                                            read_data, byref(dataId))
	assert(ier == pycf.NC_NOERR)

	# fix topology
	#ier = pycf.nccf.nccf_fix_grid_periodic_topology(gridId)

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

def inquireDataSizes(dataId):
	dims = (c_int * ndims)()
	ier = pycf.nccf.nccf_inq_data_dims(dataId, dims)
	assert(ier == pycf.NC_NOERR)
	ntot = reduce(lambda x, y: x*y, dims[:], 1)
	return ntot, dims

def getDataAsArray(dataId):
	xtypep = c_int()
	dataPtr = POINTER(c_double)()
	fillValuePtr = c_void_p()
	ier = pycf.nccf.nccf_get_data_pointer(dataId, byref(xtypep),
                                          byref(dataPtr), byref(fillValuePtr))
	assert(ier == pycf.NC_NOERR)
	assert(xtypep.value == pycf.NC_DOUBLE)
	ntot, dims = inquireDataSizes(dataId)
	data = numpy.ctypeslib.as_array(dataPtr, shape=(ntot,))
	# return a copy
	return data.copy()

def initializeData(dataId, value):
	xtypep = c_int()
	fillValuePtr = c_void_p()

	dataPtr = POINTER(c_double)()

	ier = pycf.nccf.nccf_get_data_pointer(dataId, byref(xtypep),
                                          byref(dataPtr), byref(fillValuePtr))
	assert(ier == pycf.NC_NOERR)
	ier = pycf.nccf.nccf_get_data_pointer(dataId, byref(xtypep),
                                          byref(dataPtr), byref(fillValuePtr))
	assert(ier == pycf.NC_NOERR)

	ntot, dims = inquireDataSizes(dataId)
	data = numpy.ctypeslib.as_array(dataPtr, shape=tuple(dims))
	data[...] = value

def getGridAndData(dataId):
	xtypep = c_int()
	fillValuePtr = c_void_p()
	gridId = c_int()
	coordIds = (c_int * ndims)()
	dataPtr = POINTER(c_double)()
	latPtr = POINTER(c_double)()
	lonPtr = POINTER(c_double)()
	ier = pycf.nccf.nccf_get_data_pointer(dataId, byref(xtypep),
                                          byref(dataPtr), byref(fillValuePtr))
	assert(ier == pycf.NC_NOERR)
	ier = pycf.nccf.nccf_get_data_pointer(dataId, byref(xtypep),
                                          byref(dataPtr), byref(fillValuePtr))
	assert(ier == pycf.NC_NOERR)
	ier = pycf.nccf.nccf_inq_data_gridid(dataId, byref(gridId))
	assert(ier == pycf.NC_NOERR)
	ier = pycf.nccf.nccf_inq_grid_coordids (gridId, coordIds)
	assert(ier == pycf.NC_NOERR)
	ier = pycf.nccf.nccf_get_coord_data_pointer(coordIds[0], byref(latPtr))
	assert(ier == pycf.NC_NOERR)
	ier = pycf.nccf.nccf_get_coord_data_pointer(coordIds[1], byref(lonPtr))
	assert(ier == pycf.NC_NOERR)
	ntot, dims = inquireDataSizes(dataId)
	data = numpy.ctypeslib.as_array(dataPtr, shape=tuple(dims))
	lats = numpy.ctypeslib.as_array(latPtr, shape=tuple(dims))
	lons = numpy.ctypeslib.as_array(lonPtr, shape=tuple(dims))

	return lats, lons, data


def printInvalidDataPoints(dataId, fillValue):
	lats, lons, data = getGridAndData(dataId)
	badLats = lats[data == fillValue]
	badLons = lons[data == fillValue]
	for i in range(len(badLats)):
		print('invalid point: lat = {:.10f} lon = {:.10f}'.format(badLats[i], badLons[i]))

def plotData(dataId):
	from matplotlib import pylab
	lats, lons, data = getGridAndData(dataId)
	pylab.pcolor(lons, lats, data)
	pylab.show()


timeStats = {
	'weights': float('nan'),
	'evaluation': float('nan'),
}

srcGridId, srcDataId = createData(src_file, b"src")
dstGridId, dstDataId = createData(dst_file, b"dst")

# compute the interpolation weights
regridId = c_int()
ier = pycf.nccf.nccf_def_regrid(srcGridId, dstGridId, byref(regridId))
assert(ier == pycf.NC_NOERR)
nitermax = c_int(args.nitermax)
tolpos = c_double(args.tolpos)

tic = time.time()
ier = pycf.nccf.nccf_compute_regrid_weights(regridId,
                                            nitermax, tolpos)
toc = time.time()
assert(ier == pycf.NC_NOERR)
timeStats['weights'] = toc - tic

# get the the number of valid target points
nvalid = c_int()
ier = pycf.nccf.nccf_inq_regrid_nvalid(regridId, byref(nvalid))
assert(ier == pycf.NC_NOERR)

# store the reference data values
dstDataRef = getDataAsArray(dstDataId)

# initialize the data
initializeData(dstDataId, -2.0)

# interpolate
tic = time.time()
ier = pycf.nccf.nccf_apply_regrid(regridId, srcDataId, dstDataId)
toc = time.time()
assert(ier == pycf.NC_NOERR)
timeStats['evaluation'] = toc - tic

dstDataInterp = getDataAsArray(dstDataId)

srcNtot, srcDims = inquireDataSizes(srcDataId)
dstNtot, dstDims = inquireDataSizes(dstDataId)

# compute error
error =  numpy.sum(abs(dstDataInterp - dstDataRef)) / float(dstNtot)
print('libcf interpolation:')
print('\tsrc: {} ntot: {}'.format(srcDims[:], srcNtot))
print('\tdst: {} ntot: {}'.format(dstDims[:], dstNtot))
ninvalid = dstNtot - nvalid.value
print('\t     # invalid points: {} ({:.3f}%)'.format(ninvalid,
	                                           100*ninvalid/float(dstNtot)))

printInvalidDataPoints(dstDataId, fillValue=-2.0)


print('interpolation error: {:.3g}'.format(error))
print('time stats:')
totTime = 0.0
for k, v in timeStats.items():
	print('\t{0:<32} {1:>.3g} sec'.format(k, v))
	totTime += v
print('\t{0:<32} {1:>.3g} sec'.format('total', totTime))

plotData(dstDataId)

# clean up
destroyData(srcDataId)
destroyData(dstDataId)
