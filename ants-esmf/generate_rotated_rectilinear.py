import ants.tests.stock as stock
import copy
import iris
import ants
import cartopy.crs as ccrs
import numpy
import argparse
from ants.regrid._esmf import _LatLonExtractor


parser = argparse.ArgumentParser(description='Generate rotated grid in 2d')
parser.add_argument('--nj', type=int, dest='nj', default=101,
                    help='Latitude axis dimension')
parser.add_argument('--ni', type=int, dest='ni', default=201,
                    help='Longitude axis dimension')
parser.add_argument('--pole_lat', type=float, dest='pole_lat', default=70.0,
                    help='Pole latitude')
parser.add_argument('--pole_lon', type=float, dest='pole_lon', default=20.0,
                    help='Pole longitude')
parser.add_argument('--output', type=str, dest='output', default='rotated.nc',
                    help='Output file')
args = parser.parse_args()

cube = stock.geodetic((args.nj, args.ni), with_bounds=True)
crs = iris.coord_systems.RotatedGeogCS(grid_north_pole_latitude=args.pole_lat,
	                                   grid_north_pole_longitude=args.pole_lon)
x_coord, y_coord = ants.utils.cube.horizontal_grid(cube)

x_coord.coord_system = crs
x_coord.standard_name = 'grid_longitude'
x_coord.units = 'degrees'

y_coord.coord_system = crs
y_coord.standard_name = 'grid_latitude'
y_coord.units = 'degrees'

extractor = _LatLonExtractor(cube, staggering='center')
lats = extractor.get_latitude()
lons = extractor.get_longitude()

cube.data = (numpy.sin(2*lons*numpy.pi/180.) * numpy.cos(lats*numpy.pi/180.))**2
cube.standard_name = 'air_temperature'
cube.var_name = 'cellData'
cube.grid_mapping = 'rotated_latitude_longitude'

iris.save(cube, args.output)