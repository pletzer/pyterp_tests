from matplotlib import pylab
import argparse
import iris
import numpy

iris.FUTURE.netcdf_promote = True

parser = argparse.ArgumentParser(description='Plot grid')
parser.add_argument('--src_file', type=str, dest='src_file', default='src.nc',
                    help='Source data file name')
parser.add_argument('--dst_file', type=str, dest='dst_file', default='dst.nc',
                    help='Destination data file name')

args = parser.parse_args()

pylab.figure()

def plotCellAreas(cube):
	coords = cube.coords()
	yy = coords[1].points
	xx = coords[0].points
	dyy = yy[1:, :] - yy[:-1, :]
	dxx = xx[:, 1:] - xx[:, :-1]
	dyyMid = 0.5*(dyy[:, :-1] + dyy[:, 1:])
	dxxMid = 0.5*(dxx[:-1, :] + dxx[1:, :])
	areas = numpy.zeros(yy.shape, yy.dtype)
	areas[:-1, :-1] = dyyMid * dxxMid
	negativeAreas = numpy.zeros(yy.shape, yy.dtype)
	nj, ni = yy.shape
	negativeAreas[numpy.where(areas < 0)] = 1.0
	pylab.pcolor(xx, yy, negativeAreas, cmap='bone_r')

def plotData(cube):
	coords = cube.coords()
	xx = coords[0].points
	yy = coords[1].points
	pylab.pcolor(xx, yy, cube.data)	

def plotGrid(cube, lineType='k--'):
	coords = cube.coords()
	xx = coords[0].points
	yy = coords[1].points
	nj, ni = yy.shape
	for j in range(0, nj):
		y = yy[j, :]
		x = xx[j, :]
		pylab.plot(x, y, lineType)
	for i in range(0, ni):
		y = yy[:, i]
		x = xx[:, i]
		pylab.plot(x, y, lineType)


srcCubes = iris.load(args.src_file)
srcPointCube = None
for cb in srcCubes:
    if cb.var_name == 'pointData':
        srcPointCube = cb
    if cb.var_name == 'cellData':
    	srcCellCube = cb

dstCubes = iris.load(args.dst_file)
dstPointCube = None
dstCellCube = None
for cb in dstCubes:
    if cb.var_name == 'pointData':
        dstPointCube = cb
    if cb.var_name == 'cellCube':
    	dstCellCube = cb

plotCellAreas(srcPointCube)
plotData(srcPointCube)
plotGrid(srcPointCube, 'g-')
plotGrid(dstPointCube, 'r-')

pylab.show()
