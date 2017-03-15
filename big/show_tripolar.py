import vtk
import netCDF4
import numpy

nc = netCDF4.Dataset('coords_CF_ORCA12_GO6-2.nc')
lats = nc.variables['latt'][:]
lons = nc.variables['lont'][:]

n0, n1 = lats.shape

nlines = 10

def spherePointsFromLatLons(lats, lons, radius=1.0):
	xyz = numpy.zeros((len(lats), 3), numpy.float64)
	xyz[:, 0] = radius * numpy.cos(lats * numpy.pi/180.) * numpy.cos(lons * numpy.pi/180.)
	xyz[:, 1] = radius * numpy.cos(lats * numpy.pi/180.) * numpy.sin(lons * numpy.pi/180.)
	xyz[:, 2] = radius * numpy.sin(lats * numpy.pi/180.)
	return xyz

# create nlines coordinate curves
pipeline = {'stuff': [],
            'actors': []}
for j in range(0, n0, n0//nlines):

	lat = lats[j, :]
	lon = lons[j, :]
	xyz = spherePointsFromLatLons(lat, lon)

	npts = len(lat)

	vxyz = vtk.vtkDoubleArray()
	vxyz.SetNumberOfComponents(3)
	vxyz.SetNumberOfTuples(npts)
	vxyz.Allocate(npts)
	vxyz.SetVoidArray(xyz, npts*3, 1)

	pts = vtk.vtkPoints()
	#pts.SetNumberOfPoints(len(lat))
	pts.SetDataTypeToDouble()
	pts.SetData(vxyz)

	line = vtk.vtkPolyLineSource()
	line.SetPoints(pts)

	mapper = vtk.vtkPolyDataMapper()
	mapper.SetInputConnection(line.GetOutputPort())

	pipeline['stuff'] += [xyz, vxyz, pts, line, mapper]

	actor = vtk.vtkActor()
	actor.SetMapper(mapper)
	#actor.GetProperty().SetColor(1., 0., 0.)

	pipeline['actors'].append(actor)

# rendering
ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren)
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)

# Add the actors to the renderer, set the background and size
for a in pipeline['actors']:
	ren.AddActor(a)

ren.SetBackground(0.5, 0.5, 0.5)
renWin.SetSize(400, 400)

# Interact with the data.
iren.Initialize()
renWin.Render()
iren.Start()