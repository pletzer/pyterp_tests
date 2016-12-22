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

def generateCoordsAndData(nj, ni,
                          rhoMin=0.0, rhoMax=1.0,
                          theMin=0.0, theMax=2*math.pi):
    """
    Generate coordinate and point/cell data
    @param nj number of y
    @param ni number of x
    @param rhoMin min radius
    @param rhoMax max radius
    @param theMin min angle
    @param theMax max angle 
    @return CubeList object containing point and cell data
    """

    # generate the axes
    rho = numpy.linspace(rhoMin, rhoMax, nj)
    the = numpy.linspace(theMin, theMax, ni)

    # set the curvilinear coords and field
    yy, xx = grid_mapper.createCoords(rho, the, radius=rhoMax,)
    print('min/max x: {} {}'.format(xx.min(), xx.max()))
    print('min/max y: {} {}'.format(yy.min(), yy.max()))
    pointData = grid_mapper.createPointData(yy, xx)
    yyCells, xxCells, cellData = grid_mapper.createCellData(yy, xx)

    pointCube = iris.cube.Cube(pointData, var_name='pointData', 
                               standard_name='air_temperature', cell_methods=None)
    yyCoord = iris.coords.AuxCoord(yy, var_name='yy')
    xxCoord = iris.coords.AuxCoord(xx, var_name='xx')
    pointCube.add_aux_coord(yyCoord, data_dims=(0, 1))
    pointCube.add_aux_coord(xxCoord, data_dims=(0, 1))
    
    cellCube = iris.cube.Cube(cellData, var_name='cellData', standard_name='air_temperature')
    yyBounds, yyMid = createBoundsArray(yy)
    xxBounds, xxMid = createBoundsArray(xx)
    cellAuxYY = iris.coords.AuxCoord(yyMid, var_name='yyMid', bounds=yyBounds)
    cellAuxXX = iris.coords.AuxCoord(xxMid, var_name='xxMid', bounds=xxBounds)
    cellCube.add_aux_coord(cellAuxYY, data_dims=(0, 1))
    cellCube.add_aux_coord(cellAuxXX, data_dims=(0, 1))

    return iris.cube.CubeList([pointCube, cellCube])

def generateCoordsAndDataRectilinear(nj, ni, ymin, ymax, xmin, xmax):
    """
    Generate coordinate and point/cell data
    @param nj number of y
    @param ni number of x
    @param ymin min y 
    @param ymax max y
    @param xmin min x
    @param xmax max x
    @return CubeList object containing point and cell data
    """

    # generate the axes
    y = numpy.linspace(ymin, ymax, nj)
    x = numpy.linspace(xmin, xmax, ni)

    # set the curvilinear coords and field
    yy = numpy.outer(y, numpy.ones((ni,), numpy.float64))
    xx = numpy.outer(numpy.ones((nj,), numpy.float64), x)
    print('min/max x: {} {}'.format(xx.min(), xx.max()))
    print('min/max y: {} {}'.format(yy.min(), yy.max()))
    pointData = grid_mapper.createPointData(yy, xx)
    yyCells, xxCells, cellData = grid_mapper.createCellData(yy, xx)

    pointCube = iris.cube.Cube(pointData, var_name='pointData', 
                               standard_name='air_temperature', cell_methods=None)
    yyCoord = iris.coords.AuxCoord(yy, var_name='yy')
    xxCoord = iris.coords.AuxCoord(xx, var_name='xx')
    pointCube.add_aux_coord(yyCoord, data_dims=(0, 1))
    pointCube.add_aux_coord(xxCoord, data_dims=(0, 1))
    
    cellCube = iris.cube.Cube(cellData, var_name='cellData', standard_name='air_temperature')
    yyBounds, yyMid = createBoundsArray(yy)
    xxBounds, xxMid = createBoundsArray(xx)
    cellAuxYY = iris.coords.AuxCoord(yyMid, var_name='yyMid', bounds=yyBounds)
    cellAuxXX = iris.coords.AuxCoord(xxMid, var_name='xxMid', bounds=xxBounds)
    cellCube.add_aux_coord(cellAuxYY, data_dims=(0, 1))
    cellCube.add_aux_coord(cellAuxXX, data_dims=(0, 1))

    return iris.cube.CubeList([pointCube, cellCube])

parser = argparse.ArgumentParser(description='Generate uniform data in 2d')
parser.add_argument('--src_nj', type=int, dest='src_nj', default=101, 
                    help='Source radial dimension')
parser.add_argument('--src_ni', type=int, dest='src_ni', default=201, 
                    help='Source poloidal dimension')
parser.add_argument('--dst_nj', type=int, dest='dst_nj', default=21, 
                    help='Destination radial dimension')
parser.add_argument('--dst_ni', type=int, dest='dst_ni', default=41, 
                    help='Destination poloidal dimension')
parser.add_argument('--radius', type=float, dest='radius', default=1.0, 
                    help='Radius')
parser.add_argument('--src_file', type=str, dest='src_file', default='src.nc',
                    help='Source data file name')
parser.add_argument('--dst_file', type=str, dest='dst_file', default='dst.nc',
                    help='Destination data file name')
parser.add_argument('--dst_ymin', type=float, dest='dst_ymin', default=-1/math.sqrt(2.),
                    help='Min y on destination grid')
parser.add_argument('--dst_ymax', type=float, dest='dst_ymax', default=+1/math.sqrt(2.),
                    help='Max y value on destination grid')
parser.add_argument('--dst_xmin', type=float, dest='dst_xmin', default=-1/math.sqrt(2.),
                    help='Min x on destination grid')
parser.add_argument('--dst_xmax', type=float, dest='dst_xmax', default=1/math.sqrt(2.),
                    help='Max x on destination grid')

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
srcCubes = generateCoordsAndData(args.src_nj, args.src_ni,
                                 rhoMin=0.0, rhoMax=args.radius,
                                 theMin=0.0, theMax=2*math.pi)
iris.save(srcCubes, args.src_file)
dstCubes = generateCoordsAndDataRectilinear(args.dst_nj, args.dst_ni,
                                 ymin=args.dst_ymin, ymax=args.dst_ymax,
                                 xmin=args.dst_xmin, xmax=args.dst_xmax)
iris.save(dstCubes, args.dst_file)
