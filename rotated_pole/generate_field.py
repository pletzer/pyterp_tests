import argparse
import numpy
import iris
import sys
import math
import grid_mapper

iris.FUTURE.netcdf_no_unlimited = True

parser = argparse.ArgumentParser(description='Generate uniform data in 2d')
parser.add_argument('--src_nj', type=int, dest='src_nj', default=41, 
                    help='Source latitude axis dimension')
parser.add_argument('--src_ni', type=int, dest='src_ni', default=81, 
                    help='Source longitude axis dimension')
parser.add_argument('--dst_nj', type=int, dest='dst_nj', default=121, 
                    help='Destination latitude axis dimension')
parser.add_argument('--dst_ni', type=int, dest='dst_ni', default=241, 
                    help='Destination longitude axis dimension')
parser.add_argument('--delta_theta', type=float, dest='delta_theta', default=30.0, 
                    help='Latitude of displaced pole')
parser.add_argument('--delta_lambda', type=float, dest='delta_lambda', default=00.0, 
                    help='Longitude of displaced pole')
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


# generate the axes
srcLatsPrime = numpy.linspace(latMin, latMax, args.src_nj)
srcLonsPrime = numpy.linspace(lonMin, lonMax, args.src_ni)
dstLatsPrime = numpy.linspace(latMin, latMax, args.dst_nj)
dstLonsPrime = numpy.linspace(lonMin, lonMax, args.dst_ni)

# set the curvilinear coords and field
srcLats, srcLons, srcData = grid_mapper.createCoordAndData(srcLatsPrime, srcLonsPrime, 
    delta_theta=args.delta_theta, delta_lambda=args.delta_lambda)
srcLatMin, srcLatMax = min(srcLats.flat), max(srcLats.flat)
srcLonMin, srcLonMax = min(srcLons.flat), max(srcLons.flat)
print('src lat: min = {} max = {}'.format(srcLatMin, srcLatMax))
print('src lon: min = {} max = {}'.format(srcLonMin, srcLonMax))

# target grid is regular lat-lon
dstLats, dstLons, dstData = grid_mapper.createCoordAndData(dstLatsPrime, dstLonsPrime,
    delta_theta=0.0, delta_lambda=0.0)
dstLatMin, dstLatMax = min(dstLats.flat), max(dstLats.flat)
dstLonMin, dstLonMax = min(dstLons.flat), max(dstLons.flat)
print('dst lat: min = {} max = {}'.format(dstLatMin, dstLatMax))
print('dst lon: min = {} max = {}'.format(dstLonMin, dstLonMax))

srcLatCoord = iris.coords.AuxCoord(srcLats, standard_name='latitude', units='degrees_north')
srcLonCoord = iris.coords.AuxCoord(srcLons, standard_name='longitude', units='degrees_east')
srcCube = iris.cube.Cube(srcData, standard_name='air_temperature', cell_methods=None)
srcCube.add_aux_coord(srcLatCoord, data_dims=(0, 1))
srcCube.add_aux_coord(srcLonCoord, data_dims=(0, 1))

dstLatCoord = iris.coords.AuxCoord(dstLats, standard_name='latitude', units='degrees_north')
dstLonCoord = iris.coords.AuxCoord(dstLons, standard_name='longitude', units='degrees_east')
dstCube = iris.cube.Cube(dstData, standard_name='air_temperature', cell_methods=None)
dstCube.add_aux_coord(dstLatCoord, data_dims=(0, 1))
dstCube.add_aux_coord(dstLonCoord, data_dims=(0, 1))

# save the result
iris.save(srcCube, args.src_file)
iris.save(dstCube, args.dst_file)
