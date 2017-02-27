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
    cubes = iris.load(filename)
    cubePoint, cubeCell = None, None
    # find the point and cell cubes
    for cb in cubes:
        if cb.var_name == 'pointData':
            cubePoint = cb
        if cb.var_name == 'cellData':
            cubeCell = cb
    coordsPoint = cubePoint.coords()
    latsPoint = coordsPoint[0].points
    lonsPoint = coordsPoint[1].points
    data = cubeCell.data

    return latsPoint, lonsPoint, data
    

timeStats = {
    'weights': float('nan'),
    'evaluation': float('nan'),
}

srcLatsCoords, srcLonsCoords, srcData = createData(src_file, b"src")
dstLatsCoords, dstLonsCoords, dstData = createData(dst_file, b"dst")

# save the reference (exact) field data
dstDataRef = dstData.copy()

# compute the interpolation weights
tic = time.time()
interp = sigrid.conserveInterp2D.ConserveInterp2D()
interp.setDstGrid(dstLatsCoords, dstLonsCoords)
periodicity = (False, True) # NEED TO CHECK PERIODICITY
interp.setSrcGrid(periodicity, srcLatsCoords, srcLonsCoords)
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
print('\tsrc: {} ntot: {}'.format(srcLatsCoords.shape, srcNtot))
print('\tdst: {} ntot: {}'.format(dstLatsCoords.shape, dstNtot))
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

# plot
if args.plot:
    from matplotlib import pylab
    latsCell = 0.25 * (dstLatsCoords[0:-1, 0:-1] + dstLatsCoords[1:, 0:-1] + dstLatsCoords[1:, 1:] + dstLatsCoords[0:-1, 1:])
    lonsCell = 0.25 * (dstLonsCoords[0:-1, 0:-1] + dstLonsCoords[1:, 0:-1] + dstLonsCoords[1:, 1:] + dstLonsCoords[0:-1, 1:])
    pylab.pcolor(lonsCell, latsCell, dstData)
    pylab.show()
