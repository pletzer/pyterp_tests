from __future__ import print_function
import pycf
import iris
import numpy
import sys
from ctypes import byref, c_int, c_double, c_float, POINTER, c_char_p, c_void_p
import argparse
from functools import reduce
import time

parser = argparse.ArgumentParser(description='Interpolate using libcf')
parser.add_argument('--src_field', type=str, dest='src_field', default='pointData',
                    help='Name of the source field')
parser.add_argument('--src_file', type=str, dest='src_file', default='src.nc',
                    help='Source data file name')
parser.add_argument('--dst_file', type=str, dest='dst_file', default='dst.nc',
                    help='Destination data file name')
parser.add_argument('--tolpos', type=float, dest='tolpos', default=1.e-6,
                    help='Tolerance in target space')
parser.add_argument('--nitermax', type=int, dest='nitermax', default=1000,
                    help='Max number of iterations')
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
src_field = args.src_field.encode('UTF-8') # python3
ndims = 2

def createData(filename, prefix, fieldname):
    # use iris to read in the data
    # then pass the array to create libcf objects
    cubes = iris.load(filename)
    cube = None
    for cb in cubes:
        if cb.var_name == fieldname:
            cube = cb
    coords = cube.coords()
    lats = coords[0].points
    lons = coords[1].points
    
    # create coordinates
    save = 1 # copy and save
    latId, lonId = c_int(), c_int()
    dims = (c_int * 2)(lats.shape[0], lats.shape[1])
    dimNames = (c_char_p * 2)("y", "x")
    ier = pycf.nccf.nccf_def_lat_coord(ndims, dims, dimNames, lats.ctypes.data_as(POINTER(c_double)), save, byref(latId))
    assert(ier == pycf.NC_NOERR)
    ier = pycf.nccf.nccf_def_lon_coord(ndims, dims, dimNames, lons.ctypes.data_as(POINTER(c_double)), save, byref(lonId))
    assert(ier == pycf.NC_NOERR)

    # create the grid
    gridId = c_int()
    coordIds = (c_int * 2)(latId.value, lonId.value)
    gridname = prefix + 'grid'
    ier = pycf.nccf.nccf_def_grid(coordIds, gridname, byref(gridId))
    assert(ier == pycf.NC_NOERR)

    # create the data
    dataId = c_int()
    dataname = fieldname
    ier = pycf.nccf.nccf_def_data(gridId, dataname, 'temperature', 'K', None, byref(dataId))
    assert(ier == pycf.NC_NOERR)
    ier = pycf.nccf.nccf_set_data_double(dataId, cube.data.ctypes.data_as(POINTER(c_double)))
    assert(ier == pycf.NC_NOERR)

    # get a pointer to the array
    dataPtr = POINTER(c_double)()
    xtype = c_int()
    fillValuePtr = c_void_p()
    ier = pycf.nccf.nccf_get_data_pointer(dataId, byref(xtype),
                                          byref(dataPtr), byref(fillValuePtr))
    assert(ier == pycf.NC_NOERR)

    # create a numpy array from that pointer
    array = numpy.ctypeslib.as_array(dataPtr, shape=cube.data.shape)


    return {'gridId': gridId, 'dataId': dataId, 'dataArray': array, 'lats': lats, 'lons': lons}


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

def printInvalidDataPoints(lats, lons, data, fillValue):
    badLats = lats[data == fillValue]
    badLons = lons[data == fillValue]
    for i in range(len(badLats)):
        print('invalid point: lat = {:.10f} lon = {:.10f}'.format(badLats[i], badLons[i]))

def plotData(lats, lons, data):
    from matplotlib import pylab
    pylab.pcolor(lons, lats, data)
    pylab.show()


timeStats = {
    'weights': float('nan'),
    'evaluation': float('nan'),
}

src = createData(src_file, b"src", args.src_field)
dst = createData(dst_file, b"dst", args.src_field)

# compute the interpolation weights
regridId = c_int()
ier = pycf.nccf.nccf_def_regrid(src['gridId'], dst['gridId'], byref(regridId))
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
dstDataRef = dst['dataArray'].copy()

# initialize the data
dst['dataArray'][...] = -2.0

# interpolate
tic = time.time()
ier = pycf.nccf.nccf_apply_regrid(regridId, src['dataId'], dst['dataId'])
toc = time.time()
assert(ier == pycf.NC_NOERR)
timeStats['evaluation'] = toc - tic

srcDims = src['dataArray'].shape
srcNtot = srcDims[0] * srcDims[1]
dstDims = dst['dataArray'].shape
dstNtot = dstDims[0] * dstDims[1]

# compute error
error =  numpy.sum(abs(dst['dataArray'] - dstDataRef)) / float(dstNtot)
print('libcf interpolation:')
print('\tsrc: {} ntot: {}'.format(srcDims[:], srcNtot))
print('\tdst: {} ntot: {}'.format(dstDims[:], dstNtot))
ninvalid = dstNtot - nvalid.value
print('\t     # invalid points: {} ({:.3f}%)'.format(ninvalid,
                                               100*ninvalid/float(dstNtot)))

printInvalidDataPoints(dst['lats'], dst['lons'], dst['dataArray'], fillValue=-2.0)


print('interpolation error: {:.3g}'.format(error))
print('time stats:')
totTime = 0.0
for k, v in timeStats.items():
    print('\t{0:<32} {1:>.3g} sec'.format(k, v))
    totTime += v
print('\t{0:<32} {1:>.3g} sec'.format('total', totTime))

if args.plot:
    plotData(dst['lats'], dst['lons'], dst['dataArray'])

# clean up
destroyData(src['dataId'])
destroyData(dst['dataId'])
