import vtk
import netCDF4
import numpy

nc = netCDF4.Dataset('coords_CF_ORCA12_GO6-2.nc')
lats = nc.variables['latt'][:]
lons = nc.variables['lont'][:]

n0, n1 = lats.shape

nlines = 40

def spherePointsFromLatLons(lats, lons, radius=1.0):
	xyz = numpy.zeros((len(lats), 3), numpy.float64)
	xyz[:, 0] = radius * numpy.cos(lats * numpy.pi/180.) * numpy.cos(lons * numpy.pi/180.)
	xyz[:, 1] = radius * numpy.cos(lats * numpy.pi/180.) * numpy.sin(lons * numpy.pi/180.)
	xyz[:, 2] = radius * numpy.sin(lats * numpy.pi/180.)
	return xyz

# create nlines coordinate curves
pipeline = {'stuff': [],
            'actors': []}

# earth
earth = vtk.vtkEarthSource()
earth.SetRadius(0.99)
earth.OutlineOn()
tubes = vtk.vtkTubeFilter()
tubes.SetInputConnection(earth.GetOutputPort())
tubes.SetRadius(0.01)
tubes.SetNumberOfSides(5)
emapper = vtk.vtkPolyDataMapper()
emapper.SetInputConnection(tubes.GetOutputPort())
eactor = vtk.vtkActor()
eactor.GetProperty().SetColor(0.1, 0.1, 0.1)
eactor.SetMapper(emapper)

sphere = vtk.vtkSphereSource()
sphere.SetRadius(0.98)
sphere.SetThetaResolution(257)
sphere.SetPhiResolution(129)
smapper = vtk.vtkPolyDataMapper()
smapper.SetInputConnection(sphere.GetOutputPort())
sactor = vtk.vtkActor()
sactor.GetProperty().SetColor(0.6, 0.6, 0.6)
sactor.SetMapper(smapper)

pipeline['actors'] += [eactor, sactor]

def addPipeline(xyz, pipeline, color=(0., 0., 0.)):
	npts = xyz.shape[0]
	vxyz = vtk.vtkDoubleArray()
	vxyz.SetNumberOfComponents(3)
	vxyz.SetNumberOfTuples(npts)
	vxyz.SetVoidArray(xyz, npts*3, 1)

	pts = vtk.vtkPoints()
	pts.SetData(vxyz)

	line = vtk.vtkPolyLineSource()
	line.SetPoints(pts)

	mapper = vtk.vtkPolyDataMapper()
	mapper.SetInputConnection(line.GetOutputPort())

	pipeline['stuff'] += [xyz, vxyz, pts, line, mapper]

	actor = vtk.vtkActor()
	actor.SetMapper(mapper)
	actor.GetProperty().SetColor(color)

	pipeline['actors'].append(actor)

for j in range(0, n0, n0//nlines):
	lat = lats[j, :]
	lon = lons[j, :]
	xyz = spherePointsFromLatLons(lat, lon)
	addPipeline(xyz, pipeline, color=(1., 1., 0.))

for i in range(0, n1, n1//nlines):
	lat = lats[:, i]
	lon = lons[:, i]
	xyz = spherePointsFromLatLons(lat, lon)
	addPipeline(xyz, pipeline, color=(0., 1., 1.))

# rendering
ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren)
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)

# Add the actors to the renderer, set the background and size
for a in pipeline['actors']:
	ren.AddActor(a)

ren.SetBackground(0.1*135./255., 0.1*206./255., 0.3*235./255.)
renWin.SetSize(400, 400)

# Interact with the data.
iren.Initialize()
renWin.Render()
iren.Start()