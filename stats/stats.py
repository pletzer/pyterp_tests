from __future__ import print_function
import ESMF
import iris
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

parser = argparse.ArgumentParser(description='Conservatively interpolate using ESMF')
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

    coordsPoint = cubePoint.coords()
    latsPoint = coordsPoint[0].points
    lonsPoint = coordsPoint[1].points
    
    # create the ESMF grid object
    cellDims = numpy.array([latsPoint.shape[0] - 1, latsPoint.shape[1] - 1])
    grid = ESMF.Grid(max_index=cellDims, coord_sys=ESMF.api.constants.CoordSys.SPH_DEG) #, num_peri_dims=1, periodic_dim=1)
    grid.add_coords(staggerloc=ESMF.StaggerLoc.CORNER, coord_dim=LAT_INDEX)
    grid.add_coords(staggerloc=ESMF.StaggerLoc.CORNER, coord_dim=LON_INDEX)

    coordLatsPoint = grid.get_coords(coord_dim=LAT_INDEX, staggerloc=ESMF.StaggerLoc.CORNER)
    coordLonsPoint = grid.get_coords(coord_dim=LON_INDEX, staggerloc=ESMF.StaggerLoc.CORNER)

    # get the local start/end index sets and set the point coordinates
    iBeg0 = grid.lower_bounds[ESMF.StaggerLoc.CORNER][LON_INDEX]
    iEnd0 = grid.upper_bounds[ESMF.StaggerLoc.CORNER][LON_INDEX]
    iBeg1 = grid.lower_bounds[ESMF.StaggerLoc.CORNER][LAT_INDEX]
    iEnd1 = grid.upper_bounds[ESMF.StaggerLoc.CORNER][LAT_INDEX]
    coordLatsPoint[...] = latsPoint[iBeg0:iEnd0, iBeg1:iEnd1]
    coordLonsPoint[...] = lonsPoint[iBeg0:iEnd0, iBeg1:iEnd1]

    # local sizes
    nodeDims = (iEnd0 - iBeg0, iEnd1 - iBeg1)

    # create and set the field
    field = ESMF.Field(grid, staggerloc=ESMF.StaggerLoc.CENTER)
    field.data[...] = cubeCell.data[iBeg0:iEnd0, iBeg1:iEnd1]

    return grid, field, nodeDims

timeStats = {
    'weights': float('nan'),
    'evaluation': float('nan'),
}

srcGrid, srcData, srcNodeDims = createData(src_file, b"src")
srcGridSq, srcDataSq, srcNodeDimsSq = createData(src_file, b"src")
srcDataSq.data[...] *= srcDataSq.data[...] # take square
dstGrid, dstData, dstNodeDims = createData(dst_file, b"dst")
dstGridSq, dstDataSq, dstNodeDimsSq = createData(dst_file, b"dst")

# save the reference (exact) field data
dstDataRef = dstData.data.copy()
dstData.data[...] = -1
dstDataSq.data[...] = -1

# compute the interpolation weights
tic = time.time()
regrid = ESMF.Regrid(srcfield=srcData, dstfield=dstData,
                     regrid_method=ESMF.api.constants.RegridMethod.CONSERVE,
                     unmapped_action=ESMF.api.constants.UnmappedAction.IGNORE)
timeStats['weights'] = time.time() - tic

# interpolate
tic = time.time()
regrid(srcData, dstData)
regrid(srcDataSq, dstDataSq)
timeStats['evaluation'] = time.time() - tic

# compute compute statistics
sigma = numpy.sqrt(dstDataSq.data[...] - dstData.data[...]**2)


# plot
if args.plot and nprocs == 1:
    xPointSrc = srcGrid.get_coords(coord_dim=LON_INDEX, staggerloc=ESMF.StaggerLoc.CORNER)
    yPointSrc = srcGrid.get_coords(coord_dim=LAT_INDEX, staggerloc=ESMF.StaggerLoc.CORNER)
    xxCellSrc = 0.25 * (xPointSrc[0:-1, 0:-1] + xPointSrc[1:, 0:-1] + xPointSrc[1:, 1:] + xPointSrc[0:-1, 1:])
    yyCellSrc = 0.25 * (yPointSrc[0:-1, 0:-1] + yPointSrc[1:, 0:-1] + yPointSrc[1:, 1:] + yPointSrc[0:-1, 1:])

    xPoint = dstGrid.get_coords(coord_dim=LON_INDEX, staggerloc=ESMF.StaggerLoc.CORNER)
    yPoint = dstGrid.get_coords(coord_dim=LAT_INDEX, staggerloc=ESMF.StaggerLoc.CORNER)
    xxCell = 0.25 * (xPoint[0:-1, 0:-1] + xPoint[1:, 0:-1] + xPoint[1:, 1:] + xPoint[0:-1, 1:])
    yyCell = 0.25 * (yPoint[0:-1, 0:-1] + yPoint[1:, 0:-1] + yPoint[1:, 1:] + yPoint[0:-1, 1:])

    from matplotlib import pylab
    pylab.figure(0)
    p0 = pylab.pcolor(xxCellSrc, yyCellSrc, srcData.data, vmin=-1.0, vmax=1.0)
    pylab.colorbar(p0)
    pylab.title('src field')
    pylab.figure(1)
    p1 = pylab.pcolor(xxCell, yyCell, dstData.data, vmin=-1.0, vmax=1.0)
    pylab.colorbar(p1)
    pylab.title('dst field')
    pylab.figure(2)
    p2 = pylab.pcolor(xxCell, yyCell, sigma) #, vmin=-1.0, vmax=1.0)
    pylab.colorbar(p2)
    pylab.title('dst std')
    pylab.show()
