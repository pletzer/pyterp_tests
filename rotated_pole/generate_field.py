import argparse
import numpy
import iris
import sys
import math

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
parser.add_argument('--theta_pole', type=float, dest='theta_pole', default=0.0, 
                    help='Latitude of displaced pole')
parser.add_argument('--lambda_pole', type=float, dest='lambda_pole', default=0.0, 
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
srcLats = numpy.zeros((args.src_nj, args.src_ni), numpy.float64)
srcLons = numpy.zeros((args.src_nj, args.src_ni), numpy.float64)
srcData = numpy.zeros((args.src_nj, args.src_ni), numpy.float64)


dstLatsPrime = numpy.linspace(latMin, latMax, args.dst_nj)
dstLonsPrime = numpy.linspace(lonMin, lonMax, args.dst_ni)
dstLats = numpy.zeros((args.dst_nj, args.dst_ni), numpy.float64)
dstLons = numpy.zeros((args.dst_nj, args.dst_ni), numpy.float64)
dstData = numpy.zeros((args.dst_nj, args.dst_ni), numpy.float64)

# set the curvilinear coords and field
for j in range(args.src_nj):
    the = math.pi * (srcLatsPrime[j] - args.theta_pole) / 180.
    cos_the = math.cos(the)
    sin_the = math.sin(the)
    for i in range(args.src_ni):
        lam = math.pi * (srcLonsPrime[i] - args.lambda_pole) / 180.
        cos_lam = math.cos(lam)
        sin_lam = math.sin(lam)
        x = cos_the * cos_lam
        y = cos_the * sin_lam
        z = sin_the
        rho = math.sqrt(x*x + y*y)
        srcLats[j, i] = math.atan2(z, rho)
        srcLons[j, i] = math.atan2(y, x)
        # arbitrary function
        srcData[j, i] = math.sin(2*math.pi*srcLons[j, i]/180.)*numpy.cos(math.pi*srcLats[j, i]/180.)

for j in range(args.dst_nj):
    the = math.pi * (dstLatsPrime[j] - args.theta_pole) / 180.
    cos_the = math.cos(the)
    sin_the = math.sin(the)
    for i in range(args.dst_ni):
        lam = math.pi * (dstLonsPrime[i] - args.lambda_pole) / 180.
        cos_lam = math.cos(lam)
        sin_lam = math.sin(lam)
        x = cos_the * cos_lam
        y = cos_the * sin_lam
        z = sin_the
        rho = math.sqrt(x*x + y*y)
        dstLats[j, i] = math.atan2(z, rho)
        dstLons[j, i] = math.atan2(y, x)
        dstData[j, i] = math.sin(2*math.pi*dstLons[j, i]/180.)*math.cos(math.pi*dstLats[j, i]/180.)

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
