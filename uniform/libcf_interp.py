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

# get and construct the src grid
srcLatAxisId = c_int()
ier = pycf.nccf.nccf_def_axis_from_file(src_file, b"latitude", byref(srcLatAxisId))
assert(ier == pycf.NC_NOERR)
"""
srcLonAxisId = c_int()
ier = pycf.nccf.nccf_def_axis_from_file(args.src_file, "longitude", byref(srcLonAxisId))
srcAxisIds = (c_int * ndims)(srcLatAxisId.value, srcLonAxisId.value)
srcLatCoordId = c_int()
ier = pycf.nccf.nccf_def_coord_from_axes(ndims, srcAxisIds, 0, "src_lats", "latitude", "degrees_north", byref(srcLatCoordId))
srcLonCoordId = c_int()
ier = pycf.nccf.nccf_def_coord_from_axes(ndims, srcAxisIds, 1, "src_lons", "longitude", "degrees_east", byref(srcLonCoordId))
srcCoordIds = (c_int * ndims)(srcLatCoordId.value, srcLonCoordId.value)
srcGridId  = c_int()
ier = pycf.nccf.nccf_def_grid(srcCoordIds, "srcGrid", byref(srcGridId))
"""


# interpolate
regridId = c_int()
#ier = pycf.nccf.nccf_def_regrid(srcGridId.value, dstGridId.value, byref(regridId))

# compute error

# clean up