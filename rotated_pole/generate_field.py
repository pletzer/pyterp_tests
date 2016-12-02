import argparse
import numpy
import iris
import sys
import math
import grid_mapper

iris.FUTURE.netcdf_no_unlimited = True

def createBoundsArray(arr):
    n, m = arr.shape[0] - 1, arr.shape[1] - 1
    arrBounds = numpy.zeros((n, m, 4), numpy.float64)
    arrBounds[:, :, 0] = arr[:-1, :-1]
    arrBounds[:, :, 1] = arr[:-1, 1:]
    arrBounds[:, :, 2] = arr[1:, 1:]
    arrBounds[:, :, 3] = arr[1:, :-1]
    midArray = 0.25 * (arr[:-1, :-1] + arr[:-1, 1:] + arr[1:, 1:] + arr[1:, :-1])
    return arrBounds, midArray

parser = argparse.ArgumentParser(description='Generate uniform data in 2d')
parser.add_argument('--src_nj', type=int, dest='src_nj', default=401, 
                    help='Source latitude axis dimension')
parser.add_argument('--src_ni', type=int, dest='src_ni', default=801, 
                    help='Source longitude axis dimension')
parser.add_argument('--dst_nj', type=int, dest='dst_nj', default=101, 
                    help='Destination latitude axis dimension')
parser.add_argument('--dst_ni', type=int, dest='dst_ni', default=201, 
                    help='Destination longitude axis dimension')
parser.add_argument('--delta_lat', type=float, dest='delta_lat', default=30.0, 
                    help='Pole displacement in latitude')
parser.add_argument('--delta_lon', type=float, dest='delta_lon', default=20.0, 
                    help='Pole displacement in longitude')
parser.add_argument('--src_file', type=str, dest='src_file', default='src.nc',
                    help='Source data file name')
parser.add_argument('--dst_file', type=str, dest='dst_file', default='dst.nc',
                    help='Destination data file name')
parser.add_argument('--dst_lonmin', type=float, dest='dst_lonmin', default=-176.0,
                    help='Min longitude value on destination grid')
parser.add_argument('--dst_lonmax', type=float, dest='dst_lonmax', default=176.0,
                    help='Max longitude value on destination grid')
parser.add_argument('--dst_latmin', type=float, dest='dst_latmin', default=-86.0,
                    help='Min latitude value on destination grid')
parser.add_argument('--dst_latmax', type=float, dest='dst_latmax', default=86.0,
                    help='Max latitude value on destination grid')

args = parser.parse_args()

if args.src_file is '':
    print('ERROR: must provide source data file name')
    parser.print_help()
    sys.exit(1)

if args.dst_file is '':
    print('ERROR: must provide destination data file name')
    parser.print_help()
    sys.exit(1)

latMin, latMax = -90.0, +90.0
lonMin, lonMax = -180.0, 180.


# generate the axes
srcLatsPrime = numpy.linspace(latMin, latMax, args.src_nj)
srcLonsPrime = numpy.linspace(lonMin, lonMax, args.src_ni)
dstLatsPrime = numpy.linspace(args.dst_latmin, args.dst_latmax, args.dst_nj)
dstLonsPrime = numpy.linspace(args.dst_lonmin, args.dst_lonmax, args.dst_ni)

# set the curvilinear coords and field
srcLats, srcLons = grid_mapper.createCoords(srcLatsPrime, srcLonsPrime, 
                                            delta_lat=args.delta_lat,
                                            delta_lon=args.delta_lon)
srcPointData = grid_mapper.createPointData(srcLats, srcLons)
srcLatCells, srcLonCells, srcCellData = grid_mapper.createCellData(srcLats, srcLons)

srcLatMin, srcLatMax = min(srcLats.flat), max(srcLats.flat)
srcLonMin, srcLonMax = min(srcLons.flat), max(srcLons.flat)
print('src lat: min = {} max = {}'.format(srcLatMin, srcLatMax))
print('src lon: min = {} max = {}'.format(srcLonMin, srcLonMax))

# target grid is regular lat-lon
dstLats, dstLons = grid_mapper.createCoords(dstLatsPrime, dstLonsPrime,
                                            delta_lat=0.0, delta_lon=0.0)
dstPointData = grid_mapper.createPointData(dstLats, dstLons)
dstLatCells, dstLonCells, dstCellData = grid_mapper.createCellData(dstLats, dstLons)

dstLatMin, dstLatMax = min(dstLats.flat), max(dstLats.flat)
dstLonMin, dstLonMax = min(dstLons.flat), max(dstLons.flat)
print('dst lat: min = {} max = {}'.format(dstLatMin, dstLatMax))
print('dst lon: min = {} max = {}'.format(dstLonMin, dstLonMax))

srcPointCube = iris.cube.Cube(srcPointData, var_name='pointData', standard_name='air_temperature', cell_methods=None)
srcLatCoord = iris.coords.AuxCoord(srcLats, standard_name='latitude', units='degrees_north')
srcLonCoord = iris.coords.AuxCoord(srcLons, standard_name='longitude', units='degrees_east')
srcPointCube.add_aux_coord(srcLatCoord, data_dims=(0, 1))
srcPointCube.add_aux_coord(srcLonCoord, data_dims=(0, 1))
srcCellCube = iris.cube.Cube(srcCellData, var_name='cellData', standard_name='air_temperature')
srcLatBounds, srcLatMid = createBoundsArray(srcLats)
srcLonBounds, srcLonMid = createBoundsArray(srcLons)
srcCellAuxLat = iris.coords.AuxCoord(srcLatMid, var_name='srcLatMid', standard_name='latitude', bounds=srcLatBounds)
srcCellAuxLon = iris.coords.AuxCoord(srcLonMid, var_name='srcLonMid', standard_name='longitude', bounds=srcLonBounds)
srcCellCube.add_aux_coord(srcCellAuxLat, data_dims=(0, 1))
srcCellCube.add_aux_coord(srcCellAuxLon, data_dims=(0, 1))

dstPointCube = iris.cube.Cube(dstPointData, var_name='pointData', standard_name='air_temperature', cell_methods=None)
dstLatCoord = iris.coords.AuxCoord(dstLats, standard_name='latitude', units='degrees_north')
dstLonCoord = iris.coords.AuxCoord(dstLons, standard_name='longitude', units='degrees_east')
dstPointCube.add_aux_coord(dstLatCoord, data_dims=(0, 1))
dstPointCube.add_aux_coord(dstLonCoord, data_dims=(0, 1))
dstCellCube = iris.cube.Cube(dstCellData, var_name='cellData', standard_name='air_temperature')
dstLatBounds, dstLatMid = createBoundsArray(dstLats)
dstLonBounds, dstLonMid = createBoundsArray(dstLons)
dstCellAuxLat = iris.coords.AuxCoord(dstLatMid, var_name='dstLatMid', standard_name='latitude', bounds=dstLatBounds)
dstCellAuxLon = iris.coords.AuxCoord(dstLonMid, var_name='dstLonMid', standard_name='longitude', bounds=dstLonBounds)
dstCellCube.add_aux_coord(dstCellAuxLat, data_dims=(0, 1))
dstCellCube.add_aux_coord(dstCellAuxLon, data_dims=(0, 1))

# save the result
iris.save(srcPointCube, args.src_file)
iris.save(srcCellCube, args.src_file)
iris.save(dstPointCube, args.dst_file)
iris.save(dstCellCube, args.dst_file)
