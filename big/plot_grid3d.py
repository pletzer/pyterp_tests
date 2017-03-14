import vtk
import argparse
import iris
import numpy
import math

iris.FUTURE.netcdf_promote = True


parser = argparse.ArgumentParser(description='Plot grid in 3d')
parser.add_argument('--src_file', type=str, dest='src_file', default='coords_CF_ORCA12_GO6-2.nc',
                    help='Source data file name')
parser.add_argument('--dst_file', type=str, dest='dst_file', default='dst.nc',
                    help='Destination data file name')
parser.add_argument('--src_field', type=str, dest='src_field', default='ocndept',
                    help='Destination data field name')
parser.add_argument('--dst_field', type=str, dest='dst_field', default='pointData',
                    help='Destination data field name')


args = parser.parse_args()

def readCube(filename, fieldname):
	cube = None
	cubes = iris.load(filename)
	for cb in cubes:
		if cb.var_name == fieldname:
			cube = cb
	return cube

def createPipeline(cube, pngfile, color=(1.,1.,1.), radius=1.0, show_mesh_as_surface=False, nlines=10):
	n0, n1 = cube.data.shape
	numPoints = n0 * n1
	sg = vtk.vtkStructuredGrid()
	pt = vtk.vtkPoints()
	#pt.SetNumberOfPoints(numPoints)
	coords = cube.coords()
	lats = coords[0].points
	lons = coords[1].points
        steps0 = n0 // nlines
        steps1 = n1 // nlines
	k = 0
	for i1 in range(0, n1, steps1):
		for i0 in range(0, n0, steps0):
			x = radius * math.cos(lats[i0, i1] * math.pi/180.) * math.cos(lons[i0, i1] * math.pi/180.)
			y = radius * math.cos(lats[i0, i1] * math.pi/180.) * math.sin(lons[i0, i1] * math.pi/180.)
			z = radius * math.sin(lats[i0, i1] * math.pi/180.)
			pt.InsertPoint(k, x, y, z)
			k += 1

	sg = vtk.vtkStructuredGrid()
	sg.SetDimensions(1, n0, n1)
	sg.SetPoints(pt)
	mp, ac = None, None
	# show mesh as a surface
	if show_mesh_as_surface:
		mp = vtk.vtkDataSetMapper()
                mp.SetInputData(sg)
		ac = vtk.vtkActor()
		ac.SetMapper(mp)           
	# show the grid as tubes
	ed = vtk.vtkExtractEdges()
	et = vtk.vtkTubeFilter()
	em = vtk.vtkPolyDataMapper()
	ea = vtk.vtkActor()
	et.SetRadius(0.01)
	ed.SetInputData(sg)
	et.SetInputConnection(ed.GetOutputPort())
	em.SetInputConnection(et.GetOutputPort())
	ea.SetMapper(em)
	ea.GetProperty().SetColor(color)
	return [ea], em, et, ed, mp, sg


def render(actors):
	# rendering stuff
	renderer = vtk.vtkRenderer()
	renderWindow = vtk.vtkRenderWindow()
	renderWindow.AddRenderer(renderer)
	renderWindowInteractor = vtk.vtkRenderWindowInteractor()
	renderWindowInteractor.SetRenderWindow(renderWindow)
	for a in actors:
		if a is not None:
			renderer.AddActor(a)
	renderer.SetBackground(.0, .0, .0)
	renderWindow.Render()
	renderWindowInteractor.Start()
		
src_cube = readCube(args.src_file, args.src_field)
src_pipeline = createPipeline(src_cube, pngfile='src.png', color=(0.5, 0.5, 0.5), radius=0.99, show_mesh_as_surface=False)

#dst_cube = readCube(args.dst_file, args.dst_field)
#dst_pipeline = createPipeline(dst_cube, pngfile='dst.png', color=(1.0, 0.0, 0.0), radius=1.05)

render(src_pipeline[0]) # + dst_pipeline[0])


