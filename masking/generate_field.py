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

def generateCoordsAndData(nj, ni, delta_lat, delta_lon, 
                         latMin, latMax, lonMin, lonMax):
    """
    Generate coordinate and point/cell data
    @param nj number of latitudes
    @param ni number of longitudes
    @param delta_lat rotated pole shift in latitude
    @param delta_lon rotated pole shift in longitude
    @return CubeList object containing point and cell data
    """

    # generate the axes
    latsPrime = numpy.linspace(latMin, latMax, nj)
    lonsPrime = numpy.linspace(lonMin, lonMax, ni)

    # set the curvilinear coords and field
    lats, lons = grid_mapper.createCoords(latsPrime, lonsPrime, 
                                          delta_lat=delta_lat,
                                          delta_lon=delta_lon)
    print('min/max lats: {} {}'.format(lats.min(), lats.max()))
    print('min/max lons: {} {}'.format(lons.min(), lons.max()))

    pointData = grid_mapper.createPointData(lats, lons)

    # add masking 
    pointMask = numpy.ones(pointData.shape, numpy.int)
    pointMask *= ((lats + 90.0)**2 + ((lons + 180.)/2)**2 < 160.0**2)
    pointData = numpy.ma.array(pointData, mask=pointMask)

    latCells, lonCells, cellData = grid_mapper.createCellData(lats, lons)

    pointCube = iris.cube.Cube(pointData, var_name='pointData', 
                               standard_name='air_temperature', cell_methods=None)
    print(pointCube)
    print(pointCube.data)
    print(type(pointCube.data))
    latCoord = iris.coords.AuxCoord(lats, var_name='lat',
                                    standard_name='latitude', units='degrees_north')
    lonCoord = iris.coords.AuxCoord(lons, var_name='lon',
                                    standard_name='longitude', units='degrees_east')
    pointCube.add_aux_coord(latCoord, data_dims=(0, 1))
    pointCube.add_aux_coord(lonCoord, data_dims=(0, 1))
    
    cellCube = iris.cube.Cube(cellData, var_name='cellData', standard_name='air_temperature')
    latBounds, latMid = createBoundsArray(lats)
    lonBounds, lonMid = createBoundsArray(lons)
    cellAuxLat = iris.coords.AuxCoord(latMid, var_name='latMid', 
                                      standard_name='latitude', units='degrees_north', bounds=latBounds)
    cellAuxLon = iris.coords.AuxCoord(lonMid, var_name='lonMid', 
                                      standard_name='longitude', units='degrees_east', bounds=lonBounds)
    cellCube.add_aux_coord(cellAuxLat, data_dims=(0, 1))
    cellCube.add_aux_coord(cellAuxLon, data_dims=(0, 1))

    return iris.cube.CubeList([pointCube, cellCube])

parser = argparse.ArgumentParser(description='Generate uniform data in 2d')
parser.add_argument('--src_nj', type=int, dest='src_nj', default=2, 
                    help='Source latitude axis dimension')
parser.add_argument('--src_ni', type=int, dest='src_ni', default=2, 
                    help='Source longitude axis dimension')
parser.add_argument('--dst_nj', type=int, dest='dst_nj', default=11, 
                    help='Destination latitude axis dimension')
parser.add_argument('--dst_ni', type=int, dest='dst_ni', default=21, 
                    help='Destination longitude axis dimension')
parser.add_argument('--delta_lat', type=float, dest='delta_lat', default=0.0, 
                    help='Pole displacement in latitude')
parser.add_argument('--delta_lon', type=float, dest='delta_lon', default=0.0, 
                    help='Pole displacement in longitude')
parser.add_argument('--src_file', type=str, dest='src_file', default='src.nc',
                    help='Source data file name')
parser.add_argument('--dst_file', type=str, dest='dst_file', default='dst.nc',
                    help='Destination data file name')
parser.add_argument('--dst_lonmin', type=float, dest='dst_lonmin', default=-179.0,
                    help='Min longitude value on destination grid')
parser.add_argument('--dst_lonmax', type=float, dest='dst_lonmax', default=179.0,
                    help='Max longitude value on destination grid')
parser.add_argument('--dst_latmin', type=float, dest='dst_latmin', default=-89.0,
                    help='Min latitude value on destination grid')
parser.add_argument('--dst_latmax', type=float, dest='dst_latmax', default=89.0,
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

# save the result
srcCubes = generateCoordsAndData(args.src_nj, args.src_ni, args.delta_lat, args.delta_lon,
                                 latMin=-90.0, latMax=90.0,
                                 lonMin=-180., lonMax=180.0)
iris.save(srcCubes, args.src_file)
dstCubes = generateCoordsAndData(args.dst_nj, args.dst_ni, delta_lat=0.0, delta_lon=0.0,
                                 latMin=args.dst_latmin, latMax=args.dst_latmax,
                                 lonMin=args.dst_lonmin, lonMax=args.dst_lonmax)
iris.save(dstCubes, args.dst_file)
