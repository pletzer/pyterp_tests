import vtk
import argparse
import iris
import numpy
import math

iris.FUTURE.netcdf_promote = True


parser = argparse.ArgumentParser(description='Plot grid in 3d')
parser.add_argument('--src_file', type=str, dest='src_file', default='src.nc',
                    help='Source data file name')
parser.add_argument('--dst_file', type=str, dest='dst_file', default='dst.nc',
                    help='Destination data file name')

args = parser.parse_args()

def readCube(filename):
	cube = None
	cubes = iris.load(filename)
	for cb in cubes:
		if cb.var_name == 'pointData':
			cube = cb
	return cube

def createPipeline(cube):
	cube = readCube(args.src_file)
	n0, n1 = cube.data.shape
	numPoints = n0 * n1
	sg = vtk.vtkStructuredGrid()
	pt = vtk.vtkPoints()
	pt.SetNumberOfPoints(numPoints)
	coords = cube.coords()
	lats = coords[0].points
	lons = coords[1].points
	k = 0
	for i1 in range(n1):
		for i0 in range(n0):
			x = math.cos(lats[i0, i1] * math.pi/180.) * math.cos(lons[i0, i1] * math.pi/180.)
			y = math.cos(lats[i0, i1] * math.pi/180.) * math.sin(lons[i0, i1] * math.pi/180.)
			z = math.sin(lats[i0, i1] * math.pi/180.)
			pt.SetPoint(k, x, y, z)
			k += 1

	sg = vtk.vtkStructuredGrid()
	sg.SetDimensions(1, n0, n1)
	sg.SetPoints(pt)
	mp = vtk.vtkDataSetMapper()
	mp.SetInputData(sg)
	ac = vtk.vtkActor()
	ac.SetMapper(mp)
	return ac, mp, sg, pt
		
src_cube = readCube(args.src_file)
src_pipeline = createPipeline(src_cube)

# rendering stuff
renderer = vtk.vtkRenderer()
renderWindow = vtk.vtkRenderWindow()
renderWindow.AddRenderer(renderer)
renderWindowInteractor = vtk.vtkRenderWindowInteractor()
renderWindowInteractor.SetRenderWindow(renderWindow)
renderer.AddActor(src_pipeline[0])
renderer.SetBackground(.0, .0, .0)
renderWindow.Render()
renderWindowInteractor.Start()

