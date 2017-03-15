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
ncolors = 256 
lu.SetNumberOfTableValues(ncolors)
di = 1.0 / float(ncolors - 1)
beige = numpy.array((245./255., 222./255., 179./255.))
turquoise = numpy.array((64/255., 224/255., 208/255.))
lightblue = numpy.array((0., 191./255., 255./255.))
blue = numpy.array((0., 0., 255./255.))
black = numpy.array((0., 0., 0.))
for i in range(ncolors):
	x = i * di
	if x == 0:
		rgb = beige
	elif x < 0.2:
		rgb = ((0.2 - x)/(0.2 - 0.0))*turquoise + ((x - 0.0)/(0.2 - 0.0))*lightblue
	elif x < 0.5:
		rgb = ((0.5 - x)/(0.5 - 0.2))*lightblue + ((x - 0.2)/(0.5 - 0.2))*blue
	else:
		rgb = ((1.0 - x)/(1.0 - 0.5))*blue + ((x - 0.5)/(1.0 - 0.5))*black
		
	lu.SetTableValue(i, rgb[0], rgb[1], rgb[2])

#lu.SetHueRange(0.67, 0.0)
lu.SetTableRange(0., vdata.GetMaxNorm())

mp = vtk.vtkDataSetMapper()
mp.SetInputData(sg)
mp.SetLookupTable(lu)
mp.ScalarVisibilityOn()
mp.UseLookupTableScalarRangeOn()

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
	xyz = spherePointsFromLatLons(lat, lon, radius=1.01)
	addPipeline(xyz, pipeline, color=(1., 0., 1.))

for i in range(0, n1, n1//nlines):
	lat = lats[:, i]
	lon = lons[:, i]
	xyz = spherePointsFromLatLons(lat, lon, radius=1.01)
	addPipeline(xyz, pipeline, color=(1., 0., 1.))

# rendering
ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren)
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)

# Add the actors to the renderer, set the background and size
for a in pipeline['actors']:
	ren.AddActor(a)

ren.SetBackground(0.3, 0.3, 0.3)
renWin.SetSize(860, 860)

# Interact with the data.
iren.Initialize()
renWin.Render()
iren.Start()