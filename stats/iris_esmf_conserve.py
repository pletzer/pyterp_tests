from __future__ import print_function
import ESMF
import iris
import iris.experimental.regrid_conservative
import numpy
import sys
import argparse
from functools import reduce
import time
from mpi4py import MPI

# turn on logging
esmpy = ESMF.Manager(debug=True)

# rank of this processor
pe = MPI.COMM_WORLD.Get_rank()

# number of processes
nprocs = MPI.COMM_WORLD.Get_size()

LAT_INDEX, LON_INDEX = 1, 0

parser = argparse.ArgumentParser(description='Conservatively interpolate using Iris/ESMF')
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
    # then pass the array to the ESMF API
    cubes = iris.load(filename)
    cubePoint = iris.load(filename, iris.Constraint(cube_func = lambda c: c.var_name == 'pointData'))[0]
    cubeCell = iris.load(filename, iris.Constraint(cube_func = lambda c: c.var_name == 'cellData'))[0]

    return cubeCell

timeStats = {
    'weights': float('nan'),
    'evaluation': float('nan'),
}

srcCube = createData(src_file, b"src")
dstCube = createData(dst_file, b"dst")

# save the reference (exact) field data
dstDataRef = dstCube.data.copy()
dstCube.data[...] = -1

# interpolate
tic = time.time()
dstCube = iris.experimental.regrid_conservative.regrid_conservative_via_esmpy(srcCube, dstCube)
timeStats['total'] = time.time() - tic

# compute error
srcNtot = len(srcCube.data.flat)
dstNtot = len(dstCube.data.flat)
globalSumError = numpy.sum(abs(dstCube.data - dstDataRef))

globalTimeStats = {}
for k, v in timeStats.items():
    # max value
    ts = MPI.COMM_WORLD.gather(v, root=0)
    if ts is not None:
        globalTimeStats[k] = max(ts)

if pe == 0:
    error = globalSumError / float(globalDstNtot)
    print('iris/esmf interpolation:')
    print('\tsrc: ntot: {}'.format(globalSrcNtot))
    print('\tdst: ntot: {}'.format(globalDstNtot))
    print('interpolation error: {:.3g}'.format(error))

    totTime = 0.0
    print('time stats:')
    for k, v in globalTimeStats.items():
        print('\t{0:<32} {1:>.3g} sec'.format(k, v))
        totTime += v
    print('\t{0:<32} {1:>.3g} sec'.format('total', totTime))

# plot
if args.plot and nprocs == 1:
    xxCell = dstCube.coords()[0].points
    yyCell = dstCube.coords()[1].points

    from matplotlib import pylab
    pylab.pcolor(xxCell, yyCell, dstData.data, vmin=-1.0, vmax=1.0)
    pylab.show()
