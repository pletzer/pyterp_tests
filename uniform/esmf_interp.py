from __future__ import print_function
import ESMF
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
	data = iris.load_cube(filename, 'air_temperature')
	grid = ESMF.Grid(filename=filename, filetype)


def destroyData(dataId):
	pass


timeStats = {
	'index search': float('nan'),
	'evaluation': float('nan'),
}

srcGrid, srcData = createData(src_file, b"src")
dstGrid, dstData = createData(dst_file, b"dst")

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

tic = time.time()
ier = pycf.nccf.nccf_compute_regrid_weights(regridId,
                                            nitermax, tolpos)
toc = time.time()
assert(ier == pycf.NC_NOERR)
timeStats['index search'] = toc - tic

# get the the number of valid target points
nvalid = c_int()
ier = pycf.nccf.nccf_inq_regrid_nvalid(regridId, byref(nvalid))
assert(ier == pycf.NC_NOERR)

# store the reference data values
dstDataRef = getDataAsArray(dstDataId)

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
print('\t     # valid points: {}'.format(nvalid.value))
print('interpolation error: {}'.format(error))
print('time stats:')
for k, v in timeStats.items():
	print('\t{0:<32} {1:>.3g} sec'.format(k, v))
print()

# clean up
destroyData(srcDataId)
destroyData(dstDataId)
