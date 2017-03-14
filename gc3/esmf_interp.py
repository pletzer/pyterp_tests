from __future__ import print_function
import ESMF
import netCDF4
import numpy
import sys
import argparse
from functools import reduce
import time
from mpi4py import MPI
from mpl_toolkits.basemap import Basemap

# turn on logging
esmpy = ESMF.Manager(debug=True)

# rank of this processor
pe = MPI.COMM_WORLD.Get_rank()

# number of processes
nprocs = MPI.COMM_WORLD.Get_size()

LAT_INDEX, LON_INDEX = 1, 0

parser = argparse.ArgumentParser(description='Conservatively interpolate using ESMF')
parser.add_argument('--src_file', type=str, dest='src_file', default='ac926o_10d_19850121_19850130_grid_T.nc',
                    help='Source data file name')
parser.add_argument('--src_field', type=str, dest='src_field', 
                    default='sorunoff',
                    help='Source data field name')
parser.add_argument('--dst_file', type=str, dest='dst_file', default='dst.nc',
                    help='Destination data file name')
parser.add_argument('--time', type=int, dest='time', default=0,
                    help='Time index')
parser.add_argument('--level', type=int, dest='level', default=0,
                    help='Level index')
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

def createData(filename, fieldname):

    # read the netcdf file header
    nc = netCDF4.Dataset(filename)

    var = nc.variables[fieldname]
    coordNames = var.coordinates.split(' ')

    # find the name of the curvilinear lat and lon coords
    latName, lonName = '', ''
    for c in coordNames:
        coord = nc.variables[c]
        if coord.standard_name == b'latitude':
            latName = c
        if coord.standard_name == b'longitude':
            lonName = c
    lats = nc.variables[latName][:]
    lons = nc.variables[lonName][:]

    # create the ESMF grid
    cellDims = numpy.array(lats.shape, numpy.int32) - 1
    grid = ESMF.Grid(max_index=cellDims, coord_sys=ESMF.api.constants.CoordSys.SPH_DEG) #, num_peri_dims=1, periodic_dim=1)

    # create coordinates
    grid.add_coords(staggerloc=ESMF.StaggerLoc.CORNER, coord_dim=LAT_INDEX)
    grid.add_coords(staggerloc=ESMF.StaggerLoc.CORNER, coord_dim=LON_INDEX)
    
    # get the local start/end index sets and set the point coordinates
    iBeg0 = grid.lower_bounds[ESMF.StaggerLoc.CORNER][LON_INDEX]
    iEnd0 = grid.upper_bounds[ESMF.StaggerLoc.CORNER][LON_INDEX]
    iBeg1 = grid.lower_bounds[ESMF.StaggerLoc.CORNER][LAT_INDEX]
    iEnd1 = grid.upper_bounds[ESMF.StaggerLoc.CORNER][LAT_INDEX]
    
    coordLatsPoint = grid.get_coords(coord_dim=LAT_INDEX, staggerloc=ESMF.StaggerLoc.CORNER)
    coordLonsPoint = grid.get_coords(coord_dim=LON_INDEX, staggerloc=ESMF.StaggerLoc.CORNER)

    # set the ESMF coordinates
    coordLatsPoint[:] = lats
    coordLonsPoint[:] = lons

    # create and set the field, cell centred
    field = ESMF.Field(grid, staggerloc=ESMF.StaggerLoc.CORNER)

    # read the cell centred data and set the field. Note that we need to use the point dims
    data = var[:]
    if len(data.shape) == 2:
        field.data[...] = data[iBeg0:iEnd0, iBeg1:iEnd1]
    elif len(data.shape) == 3:
        field.data[...] = data[args.time, iBeg0:iEnd0, iBeg1:iEnd1]
    else:
        field.data[...] = data[args.time, args.level, iBeg0:iEnd0, iBeg1:iEnd1]

    return grid, field

timeStats = {
    'weights': float('nan'),
    'evaluation': float('nan'),
}

srcGrid, srcData = createData(src_file, args.src_field)
dstGrid, dstData = createData(dst_file, 'pointData')

# initialize the dst data
dstData.data[...] = 0

# compute the interpolation weights
tic = time.time()
regrid = ESMF.Regrid(srcfield=srcData, dstfield=dstData,
                     regrid_method=ESMF.api.constants.RegridMethod.BILINEAR,
                     unmapped_action=ESMF.api.constants.UnmappedAction.IGNORE,
                     ignore_degenerate=True)
timeStats['weights'] = time.time() - tic

# interpolate
tic = time.time()
regrid(srcData, dstData)
timeStats['evaluation'] = time.time() - tic

# plot
if args.plot and nprocs == 1:

    xPoint = dstGrid.get_coords(coord_dim=LON_INDEX, staggerloc=ESMF.StaggerLoc.CORNER)
    yPoint = dstGrid.get_coords(coord_dim=LAT_INDEX, staggerloc=ESMF.StaggerLoc.CORNER)
    from matplotlib import pylab
    from mpl_toolkits.basemap import Basemap
    mp = Basemap(resolution='l')
    mp.drawcoastlines(linewidth=0.25)
    p = pylab.pcolor(xPoint, yPoint, dstData.data)
    pylab.colorbar(p)
    pylab.show()
