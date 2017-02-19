from __future__ import print_function
import ESMF
import iris
import numpy
import sys
import argparse
from functools import reduce
import time

LAT_INDEX, LON_INDEX = 1, 0

# turn on logging
esmpy = ESMF.Manager(debug=True)

parser = argparse.ArgumentParser(description='Interpolate using ESMF')
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
    cube = cube = iris.load(filename, iris.Constraint(cube_func = lambda c: c.var_name == 'pointData'))[0]
    coords = cube.coords()
    lats = coords[0].points
    lons = coords[1].points
    
    # create the ESMF grid object

    cellDims = numpy.array([lats.shape[0] - 1, lats.shape[1] - 1])
    grid = ESMF.Grid(max_index=cellDims, coord_sys=ESMF.api.constants.CoordSys.SPH_DEG)#, num_peri_dims=1) #, periodic_dim=1)

    grid.add_coords(staggerloc=ESMF.StaggerLoc.CORNER, coord_dim=LAT_INDEX)
    grid.add_coords(staggerloc=ESMF.StaggerLoc.CORNER, coord_dim=LON_INDEX)

    coordLat = grid.get_coords(coord_dim=LAT_INDEX, staggerloc=ESMF.StaggerLoc.CORNER)
    coordLon = grid.get_coords(coord_dim=LON_INDEX, staggerloc=ESMF.StaggerLoc.CORNER)

    # get the local start/end index sets
    iBegLat = grid.lower_bounds[ESMF.StaggerLoc.CORNER][LAT_INDEX]
    iEndLat = grid.upper_bounds[ESMF.StaggerLoc.CORNER][LAT_INDEX]
    iBegLon = grid.lower_bounds[ESMF.StaggerLoc.CORNER][LON_INDEX]
    iEndLon = grid.upper_bounds[ESMF.StaggerLoc.CORNER][LON_INDEX]

    # set the coordinates
    coordLat[...] = lats[iBegLon:iEndLon, iBegLat:iEndLat]
    coordLon[...] = lons[iBegLon:iEndLon, iBegLat:iEndLat]
   
    # create field
    field = ESMF.Field(grid, name="air_temperature", 
                   staggerloc=ESMF.StaggerLoc.CORNER)
    field.data[...] = cube.data[:]

    nodeDims = (iEndLat - iBegLat, iEndLon - iBegLon)

    return grid, field, nodeDims

timeStats = {
    'weights': float('nan'),
    'evaluation': float('nan'),
}

srcGrid, srcData, srcNodeDims = createData(src_file, b"src")
dstGrid, dstData, dstNodeDims = createData(dst_file, b"dst")

# save the reference (exact) field data
dstDataRef = dstData.data.copy()
dstData.data[...] = -1

# compute the interpolation weights
tic = time.time()
regrid = ESMF.api.regrid.Regrid(srcData, dstData,
                                src_mask_values=None, dst_mask_values=None,
                                regrid_method=ESMF.api.constants.RegridMethod.BILINEAR,
                                pole_method=None,
                                regrid_pole_npoints=None, # only relevant if method is ALLAVG
                                line_type=ESMF.api.constants.LineType.GREAT_CIRCLE, # how the distance between two points is computed
                                norm_type=None, # only for conservative regridding
                                unmapped_action=ESMF.api.constants.UnmappedAction.IGNORE, 
                                ignore_degenerate=True, # produce an error if two points are degenerate and if set to False
                                src_frac_field=None, dst_frac_field=None)
timeStats['weights'] = time.time() - tic

# interpolate
tic = time.time()
regrid(srcData, dstData)
timeStats['evaluation'] = time.time() - tic

# compute error
srcNtot = len(srcData.data.flat)
dstNtot = len(dstData.data.flat)
error =  numpy.sum(abs(dstData.data - dstDataRef)) / float(dstNtot)
print('emsf interpolation:')
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
if args.plot:
    lats = dstGrid.get_coords(coord_dim=LAT_INDEX, staggerloc=ESMF.StaggerLoc.CORNER)
    lons = dstGrid.get_coords(coord_dim=LON_INDEX, staggerloc=ESMF.StaggerLoc.CORNER)
    from matplotlib import pylab
    pylab.pcolor(lons, lats, dstData.data)
    pylab.show()

# clean up
# nothing to do
