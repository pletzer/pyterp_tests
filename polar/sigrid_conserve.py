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
    cubes = iris.load(filename)
    cubePoint, cubeCell = None, None
    # find the point and cell cubes
    for cb in cubes:
        if cb.var_name == 'pointData':
            cubePoint = cb
        if cb.var_name == 'cellData':
            cubeCell = cb
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
dstData[...] = -1

# compute the interpolation weights
tic = time.time()
interp = sigrid.conserveInterp2D.ConserveInterp2D()
interp.setDstGrid(dstXCoords, dstYCoords)
interp.setSrcGrid(periodicity, srcXCoords, srcYCoords)
interp.computeWeights()
timeStats['weights'] = time.time() - tic

# interpolate
tic = time.time()
dstData = interp.apply(srcData)
timeStats['evaluation'] = time.time() - tic

# compute error
srcNtot = len(srcData.data.flat)
dstNtot = len(dstData.data.flat)
error =  numpy.sum(abs(dstData - dstDataRef)) / float(dstNtot)
print('sigrid interpolation:')
print('\tsrc: {} ntot: {}'.format(srcNodeDims, srcNtot))
print('\tdst: {} ntot: {}'.format(dstNodeDims, dstNtot))
print('interpolation error: {:.3g}'.format(error))
totTime = 0.0
print('time stats:')
for k, v in timeStats.items():
    print('\t{0:<32} {1:>.3g} sec'.format(k, v))
    totTime += v
print('\t{0:<32} {1:>.3g} sec'.format('total', totTime))

# plot
from matplotlib import pylab
pylab.pcolor(dstXCoords, dstYCoords, dstData)
pylab.show()
