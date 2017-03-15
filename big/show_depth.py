import vtk
import netCDF4
import numpy

nc = netCDF4.Dataset('coords_CF_ORCA12_GO6-2.nc')
lats = nc.variables['latt'][:]
lons = nc.variables['lont'][:]
data = nc.variables['ocndept'][:]

n0, n1 = lats.shape
npts = n0 * n1

nlines = 40

def spherePointsFromLatLons(lats, lons, radius=1.0):
	n = reduce(lambda x, y: x*y, lats.shape)
	la = lats.reshape((n,)) * numpy.pi / 180.
	lo = lons.reshape((n,)) * numpy.pi / 180.
	xyz = numpy.zeros((n, 3), numpy.float64)
	xyz[:, 0] = radius * numpy.cos(la) * numpy.cos(lo)
	xyz[:, 1] = radius * numpy.cos(la) * numpy.sin(lo)
	xyz[:, 2] = radius * numpy.sin(la)
	return xyz

# create nlines coordinate curves
pipeline = {'stuff': [],
            'actors': []}

# show depth
xyzGrid = spherePointsFromLatLons(lats, lons)

vxyzData = vtk.vtkDoubleArray()
vxyzData.SetNumberOfComponents(3)
vxyzData.SetNumberOfTuples(npts)
vxyzData.SetVoidArray(xyzGrid, npts*3, 1)

vxyz = vtk.vtkPoints()
vxyz.SetNumberOfPoints(npts)
vxyz.SetData(vxyzData)

vdata = vtk.vtkDoubleArray()
vdata.SetNumberOfComponents(1)
vdata.SetNumberOfTuples(npts)
vdata.SetVoidArray(data, npts, 1)

sg = vtk.vtkStructuredGrid()
sg.SetDimensions(n1, n0, 1) # varies faster in n1
sg.SetPoints(vxyz)
sg.GetPointData().SetScalars(vdata)

lu = vtk.vtkLookupTable()
ncolors = 16 + 1
lu.SetNumberOfTableValues(ncolors)
di = 1.0 / float(ncolors - 1)
for i in range(ncolors):
	r = max(0., 1. - 4*i*di)
	g = max(0., 1. - 2*i*di)
	b = 1. - i*di
	lu.SetTableValue(i, r, g, b)
lu.SetTableRange(0., vdata.GetMaxNorm())

mp = vtk.vtkDataSetMapper()
mp.SetInputData(sg)
mp.SetLookupTable(lu)

actor = vtk.vtkActor()
actor.SetMapper(mp)


pipeline['actors'] += [actor]

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
renWin.SetSize(860, 860)

# Interact with the data.
iren.Initialize()
renWin.Render()
iren.Start()