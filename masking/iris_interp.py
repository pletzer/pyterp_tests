from __future__ import print_function
import iris
import numpy
import sys
import argparse
import time

parser = argparse.ArgumentParser(description='Interpolate using libcf')
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

src_file = args.src_file
dst_file = args.dst_file
ndims = 2

def createData(filename, prefix):
	# use iris to read in the data
	# then pass the array to the ESMF API
	pointCube = iris.load(filename, iris.Constraint(cube_func = lambda c: c.var_name == 'pointData'))[0]
	return pointCube

timeStats = {
	'interp': float('nan'),
}

srcCube = createData(src_file, "src")
dstCube = createData(dst_file, "dst")

# save the reference (exact) field data
dstDataRef = dstCube.data.copy()

# compute the interpolation weights
scheme = iris.analysis.Linear()
tic = time.time()
interpCube = srcCube.regrid(dstCube, scheme=scheme)
timeStats['interp'] = time.time() - tic

# compute error
srcNtot = len(srcCube.data.flat)
dstNtot = len(dstCube.data.flat)
srcNodeDims = srcCube.data.shape
dstNodeDims = dstCube.data.shape
error =  numpy.sum(abs(interpCube.data - dstDataRef)) / float(dstNtot)
print('iris linear interpolation:')
print('\tsrc: {} ntot: {}'.format(srcNodeDims, srcNtot))
print('\tdst: {} ntot: {}'.format(dstNodeDims, dstNtot))
print('interpolation error: {:.3g}'.format(error))
print('time stats:')
for k, v in timeStats.items():
	print('\t{0:<32} {1:>.3g} sec'.format(k, v))
print()

# check sum
checksum = numpy.sum(dstCube.data, axis=None)
print('check sum: {:.15g}'.format(checksum))

if args.plot:
    from matplotlib import pylab
    lats = dstCube.coords()[0].points
    lons = dstCube.coords()[1].points
    data = dstCube.data
    pylab.pcolor(lons, lats, data)
    pylab.show()

# clean up
# nothing to do