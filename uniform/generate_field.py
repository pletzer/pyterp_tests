import argparse
import numpy
import iris
import sys

iris.FUTURE.netcdf_no_unlimited = True

parser = argparse.ArgumentParser(description='Generate uniform data in 2d')
parser.add_argument('--src_nj', type=int, dest='src_nj', default=101, 
                    help='Source latitude axis dimension')
parser.add_argument('--src_ni', type=int, dest='src_ni', default=201, 
                    help='Source longitude axis dimension')
parser.add_argument('--dst_nj', type=int, dest='dst_nj', default=201, 
                    help='Destination latitude axis dimension')
parser.add_argument('--dst_ni', type=int, dest='dst_ni', default=401, 
                    help='Destination longitude axis dimension')
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

latMin, latMax = -90.0, +90.0
lonMin, lonMax = 0.0, 360.0

def createPointData(lons, lats):
    # some arbitrary expression
    data = numpy.outer(numpy.cos(numpy.pi*lats/180.), numpy.sin(2*numpy.pi*lons/180.))
    return data

def createCellData(lons, lats):
    lonMid = 0.5*(lons[:-1] + lons[1:])
    latMid = 0.5*(lats[:-1] + lats[1:])
    data = numpy.outer(numpy.cos(numpy.pi*latMid/180.), numpy.sin(2*numpy.pi*lonMid/180.))
    latBounds = numpy.zeros((len(lats) - 1, 2), numpy.float64)
    latBounds[:, 0] = lats[:-1]
    latBounds[:, 1] = lats[1:]
    lonBounds = numpy.zeros((len(lons) - 1, 2), numpy.float64)
    lonBounds[:, 0] = lons[:-1]
    lonBounds[:, 1] = lons[1:]
    return {'data': data, 
            'latBounds': latBounds, 
            'lonBounds': lonBounds,
            'latMid': latMid,
            'lonMid': lonMid}

# generate the axes
srcLats = numpy.linspace(latMin, latMax, args.src_nj)
srcLons = numpy.linspace(lonMin, lonMax, args.src_ni)
dstLats = numpy.linspace(latMin, latMax, args.dst_nj)
dstLons = numpy.linspace(lonMin, lonMax, args.dst_ni)

# set the field
srcPointData = createPointData(srcLons, srcLats)
dstPointData = createPointData(dstLons, dstLats)

srcCellDict = createCellData(srcLons, srcLats)
dstCellDict = createCellData(dstLons, dstLats)

srcLatCoord = iris.coords.DimCoord(srcLats, standard_name='latitude', units='degrees_north')
srcLonCoord = iris.coords.DimCoord(srcLons, standard_name='longitude', units='degrees_east')
srcCellLatCoord = iris.coords.DimCoord(srcCellDict['latMid'], standard_name='latitude',
                                       units='degrees_north', bounds=srcCellDict['latBounds'])
srcCellLonCoord = iris.coords.DimCoord(srcCellDict['lonMid'], standard_name='longitude',
                                       units='degrees_east', bounds=srcCellDict['lonBounds'])

srcPointCube = iris.cube.Cube(srcPointData, var_name='pointData', standard_name='air_temperature')
srcPointCube.add_dim_coord(srcLatCoord, data_dim=0)
srcPointCube.add_dim_coord(srcLonCoord, data_dim=1)

srcCellCube = iris.cube.Cube(srcCellDict['data'], var_name='cellData', standard_name='air_temperature')
srcCellCube.add_dim_coord(srcCellLatCoord, data_dim=0)
srcCellCube.add_dim_coord(srcCellLonCoord, data_dim=1)

dstLatCoord = iris.coords.DimCoord(dstLats, standard_name='latitude', units='degrees_north')
dstLonCoord = iris.coords.DimCoord(dstLons, standard_name='longitude', units='degrees_east')
dstCellLatCoord = iris.coords.DimCoord(dstCellDict['latMid'], standard_name='latitude',
                                       units='degrees_north', bounds=dstCellDict['latBounds'])
dstCellLonCoord = iris.coords.DimCoord(dstCellDict['lonMid'], standard_name='longitude',
                                       units='degrees_east', bounds=dstCellDict['lonBounds'])

dstPointCube = iris.cube.Cube(dstPointData, var_name='pointData', standard_name='air_temperature')
dstPointCube.add_dim_coord(dstLatCoord, data_dim=0)
dstPointCube.add_dim_coord(dstLonCoord, data_dim=1)

dstCellCube = iris.cube.Cube(dstCellDict['data'], var_name='cellData', standard_name='air_temperature')
dstCellCube.add_dim_coord(dstCellLatCoord, data_dim=0)
dstCellCube.add_dim_coord(dstCellLonCoord, data_dim=1)

# save the result
iris.save([srcPointCube, srcCellCube], args.src_file)
iris.save([dstPointCube, dstCellCube], args.dst_file)
