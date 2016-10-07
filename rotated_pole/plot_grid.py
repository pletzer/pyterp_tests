from matplotlib import pylab
import argparse
import iris

parser = argparse.ArgumentParser(description='Plot grid')
parser.add_argument('--src_file', type=str, dest='src_file', default='src.nc',
                    help='Source data file name')
parser.add_argument('--dst_file', type=str, dest='dst_file', default='dst.nc',
                    help='Destination data file name')

args = parser.parse_args()

pylab.figure()

def plotGrid(cube, lineType='k--'):
	coords = cube.coords()
	lats = coords[0].points
	lons = coords[1].points
	nj, ni = lats.shape
	for j in range(nj):
		y = lats[j, :]
		x = lons[j, :]
		pylab.plot(x, y, lineType)
	for i in range(ni):
		y = lats[:, i]
		x = lons[:, i]
		pylab.plot(x, y, lineType)


srcCube = iris.load_cube(args.src_file, 'air_temperature')
dstCube = iris.load_cube(args.dst_file, 'air_temperature')
plotGrid(srcCube, 'g-')
plotGrid(dstCube, 'r-')

pylab.show()
