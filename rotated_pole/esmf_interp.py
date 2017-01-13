from __future__ import print_function
import ESMF
import iris
import numpy
import sys
import argparse
from functools import reduce
import time

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
    cubes = iris.load(filename)
    cube = None
    for cb in cubes:
        if cb.var_name == 'pointData':
            cube = cb
    coords = cube.coords()
    lats = coords[0].points
    lons = coords[1].points
    
    # create the ESMF grid object

    latIndex, lonIndex = 0, 1
    cellDims = numpy.array([lats.shape[0] - 1, lats.shape[1] - 1])
    grid = ESMF.Grid(max_index=cellDims) #, num_peri_dims=1, periodic_dim=1)

    grid.add_coords(staggerloc=ESMF.StaggerLoc.CORNER, coord_dim=latIndex)
    grid.add_coords(staggerloc=ESMF.StaggerLoc.CORNER, coord_dim=lonIndex)

    coordLat = grid.get_coords(coord_dim=latIndex, staggerloc=ESMF.StaggerLoc.CORNER)
    coordLon = grid.get_coords(coord_dim=lonIndex, staggerloc=ESMF.StaggerLoc.CORNER)

    # get the local start/end index sets
    iBegLat = grid.lower_bounds[ESMF.StaggerLoc.CORNER][latIndex]
    iEndLat = grid.upper_bounds[ESMF.StaggerLoc.CORNER][latIndex]
    iBegLon = grid.lower_bounds[ESMF.StaggerLoc.CORNER][lonIndex]
    iEndLon = grid.upper_bounds[ESMF.StaggerLoc.CORNER][lonIndex]

    # set the coordinates
    coordLat[...] = lats[iBegLat:iEndLat, iBegLon:iEndLon]
    coordLon[...] = lons[iBegLat:iEndLat, iBegLon:iEndLon]

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
    latIndex, lonIndex = 0, 1
    lats = dstGrid.get_coords(coord_dim=latIndex, staggerloc=ESMF.StaggerLoc.CORNER)
    lons = dstGrid.get_coords(coord_dim=lonIndex, staggerloc=ESMF.StaggerLoc.CORNER)
    from matplotlib import pylab
    pylab.pcolor(lons, lats, dstData.data)
    pylab.show()

# clean up
# nothing to do
