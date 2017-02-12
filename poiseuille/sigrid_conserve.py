from __future__ import print_function
import sigrid.conserveInterp2D
import iris
import numpy
import sys
import argparse
from functools import reduce
import time

parser = argparse.ArgumentParser(description='Conservatively interpolate using sigrid')
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
    cubePoint = iris.load(filename, iris.Constraint(cube_func = lambda c: c.var_name=='pointData'))[0]
    cubeCell = iris.load(filename, iris.Constraint(cube_func = lambda c: c.var_name=='flux'))[0]
    coordsPoint = cubePoint.coords()
    xPoint = coordsPoint[0].points
    yPoint = coordsPoint[1].points
    data = cubeCell.data

    return xPoint, yPoint, data
    

timeStats = {
    'weights': float('nan'),
    'evaluation': float('nan'),
}

srcXCoords, srcYCoords, srcData = createData(src_file, b"src")
dstXCoords, dstYCoords, dstData = createData(dst_file, b"dst")

# save the reference (exact) field data
dstDataRef = dstData.copy()

# compute the interpolation weights
tic = time.time()
interp = sigrid.conserveInterp2D.ConserveInterp2D()
interp.setDstGrid(dstXCoords, dstYCoords)
periodicity = (False, True) # NEED TO CHECK PERIODICITY
interp.setSrcGrid(periodicity, srcXCoords, srcYCoords)
interp.computeWeights()
timeStats['weights'] = time.time() - tic

# interpolate
tic = time.time()
dstData = interp.apply(srcData)
timeStats['evaluation'] = time.time() - tic

# compute error
srcNtot = len(srcData.flat)
dstNtot = len(dstData.flat)
error =  numpy.sum(abs(dstData - dstDataRef)) / float(dstNtot)
print('sigrid interpolation:')
print('\tsrc: {} ntot: {}'.format(srcXCoords.shape, srcNtot))
print('\tdst: {} ntot: {}'.format(dstXCoords.shape, dstNtot))
print('interpolation error: {:.3g}'.format(error))
totTime = 0.0
print('time stats:')
for k, v in timeStats.items():
    print('\t{0:<32} {1:>.3g} sec'.format(k, v))
    totTime += v
print('\t{0:<32} {1:>.3g} sec'.format('total', totTime))

# check sum
checksum = numpy.sum(dstData, axis=None)
print('check sum: {:.15g}'.format(checksum))

if args.plot:
    from matplotlib import pylab
    xxCell = 0.25 * (dstXCoords[0:-1, 0:-1] + dstXCoords[1:, 0:-1] + dstXCoords[1:, 1:] + dstXCoords[0:-1, 1:])
    yyCell = 0.25 * (dstYCoords[0:-1, 0:-1] + dstYCoords[1:, 0:-1] + dstYCoords[1:, 1:] + dstYCoords[0:-1, 1:])
    pylab.pcolor(xxCell, yyCell, dstData, vmin=-1.0, vmax=1.0)
    pylab.show()
