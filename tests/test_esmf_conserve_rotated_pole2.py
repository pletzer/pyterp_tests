import ESMF
import numpy
import math
from matplotlib import pylab

ESMF.Manager(debug=True)

LAT_INDEX, LON_INDEX = 1, 0

def createRotatedPoleCoords(latsPrime, lonsPrime, delta_lat, delta_lon):
    """
    Create coordinates from axes
    @param latsPrime latitude logical axis
    @param lonsPrime longitude logical axis
    @param delta_lat latitude pole displacement
    @param delta_lon longitude pole displacement
    @return curvilinear latitudes, longitudes
    """
    nj, ni = len(latsPrime), len(lonsPrime)
    shp = (ni, nj)
    lats = numpy.zeros(shp, numpy.float64)
    lons = numpy.zeros(shp, numpy.float64)

    alpha = math.pi * delta_lat / 180.
    beta = math.pi * delta_lon / 180.
    cos_alp = math.cos(alpha)
    sin_alp = math.sin(alpha)
    cos_bet = math.cos(beta)
    sin_bet = math.sin(beta)

    # http://gis.stackexchange.com/questions/10808/lon-lat-transformation
    rot_alp = numpy.array([[ cos_alp, 0., sin_alp],
                           [ 0.,      1., 0.     ],
                           [-sin_alp, 0., cos_alp]])
    rot_bet = numpy.array([[ cos_bet, sin_bet, 0.],
                           [-sin_bet, cos_bet, 0.],
                           [ 0.     , 0.,      1.]])
    transfMatrix = numpy.dot(rot_bet, rot_alp)

    xyzPrime = numpy.zeros((3,), numpy.float64)
    xyz = numpy.zeros((3,), numpy.float64)

    for j in range(nj):
        the = math.pi * latsPrime[j] / 180.
        cos_the = math.cos(the)
        sin_the = math.sin(the)
        rho = cos_the
        for i in range(ni):
            lam = math.pi * lonsPrime[i] / 180.
            cos_lam = math.cos(lam)
            sin_lam = math.sin(lam)

            xyzPrime = rho * cos_lam, rho * sin_lam, sin_the
            xyz = numpy.dot(transfMatrix, xyzPrime)

            lats[i, j] = 180. * math.asin(xyz[2]) / math.pi
            lons[i, j] = 180. * math.atan2(xyz[1], xyz[0]) / math.pi

    return lats, lons


def createGrid(nodeDims, lats2D, lons2D):
    """
    Create a structured grid from coordinates
    @param nodeDims number of latitude/longitude nodes
    @param lats2D latitude coordinates (corner points)
    @param lons2D langitude coordinates (corner points)
    @return ESMF Grid
    """
    cellDims = numpy.array(nodeDims) - 1

    grid = ESMF.Grid(cellDims,
                     num_peri_dims=0, # number of periodic dimensions
                     periodic_dim=None, # 1, 2, or 3
                     pole_dim=None, # default is 1 (1-based indexing?). Not sure if this follows Fortran ordering
                     coord_sys=ESMF.api.constants.CoordSys.SPH_DEG, # other choices are CART, SPH_RAD
                     coord_typekind=ESMF.api.constants.TypeKind.R8,
                     staggerloc=None) #ESMF.api.constants.StaggerLoc.CORNER) #not sure why this has to be set here?

    # attach nodes to the grid
    grid.add_coords(staggerloc=ESMF.StaggerLoc.CORNER, coord_dim=LAT_INDEX)
    coordLatPoint = grid.get_coords(coord_dim=LAT_INDEX, staggerloc=ESMF.StaggerLoc.CORNER)

    grid.add_coords(staggerloc=ESMF.StaggerLoc.CORNER, coord_dim=LON_INDEX)
    coordLonPoint = grid.get_coords(coord_dim=LON_INDEX, staggerloc=ESMF.StaggerLoc.CORNER)

    iBegLat = grid.lower_bounds[ESMF.StaggerLoc.CORNER][LAT_INDEX]
    iEndLat = grid.upper_bounds[ESMF.StaggerLoc.CORNER][LAT_INDEX]
    iBegLon = grid.lower_bounds[ESMF.StaggerLoc.CORNER][LON_INDEX]
    iEndLon = grid.upper_bounds[ESMF.StaggerLoc.CORNER][LON_INDEX]

    coordLatPoint[...] = lats2D[iBegLon:iEndLon, iBegLat:iEndLat]
    coordLonPoint[...] = lons2D[iBegLon:iEndLon, iBegLat:iEndLat]

    return grid

def createField(grid, name, data):
    """
    Create cell-centered structured field
    @param grid ESMF grid
    @param name field name
    @param data 2D data of size number of lat cells times number of lon cells
    @return ESMF Field
    """
    field = ESMF.Field(grid,
                       name=name,
                       typekind=ESMF.api.constants.TypeKind.R8,
                       staggerloc=ESMF.StaggerLoc.CENTER,
                       ndbounds=1) # scalar
    field.data[...] = data
    return field

def plotField(field, linetype='k-'):
    iBegLat = field.grid.lower_bounds[ESMF.StaggerLoc.CORNER][LAT_INDEX]
    iEndLat = field.grid.upper_bounds[ESMF.StaggerLoc.CORNER][LAT_INDEX]
    iBegLon = field.grid.lower_bounds[ESMF.StaggerLoc.CORNER][LON_INDEX]
    iEndLon = field.grid.upper_bounds[ESMF.StaggerLoc.CORNER][LON_INDEX]

    lats = field.grid.get_coords(LAT_INDEX, ESMF.StaggerLoc.CORNER)[iBegLon:iEndLon, iBegLat:iEndLat]
    lons = field.grid.get_coords(LON_INDEX, ESMF.StaggerLoc.CORNER)[iBegLon:iEndLon, iBegLat:iEndLat]

    # plot the grid
    nc = 20 # number of curves
    ni, nj = lats.shape
    stepj, stepi = max(1, nj//nc), max(1, ni//nc)
    for j in range(0, nj, stepj):
        pylab.plot(lons[:, j], lats[:, j], linetype)
    for i in range(0, ni, stepi):
        pylab.plot(lons[i, :], lats[i, :], linetype)

    pylab.pcolormesh(lons, lats, field.data[0,...], vmin=-180.0, vmax=180.0)
    #pylab.show()

def plotGrid3d(field, color=(0.5, 0.1, 0.1), radius=1.0):
    import vtk
    pt = vtk.vtkPoints()
    sg = vtk.vtkStructuredGrid()
    mp = vtk.vtkDataSetMapper()
    ac = vtk.vtkActor()
    # connect
    sg.SetPoints(pt)
    mp.SetInputData(sg)
    ac.SetMapper(mp)
    # set
    ac.GetProperty().SetColor(color)
    iBegLat = field.grid.lower_bounds[ESMF.StaggerLoc.CORNER][LAT_INDEX]
    iEndLat = field.grid.upper_bounds[ESMF.StaggerLoc.CORNER][LAT_INDEX]
    iBegLon = field.grid.lower_bounds[ESMF.StaggerLoc.CORNER][LON_INDEX]
    iEndLon = field.grid.upper_bounds[ESMF.StaggerLoc.CORNER][LON_INDEX]
    lats = field.grid.get_coords(LAT_INDEX, ESMF.StaggerLoc.CORNER)[iBegLon:iEndLon, iBegLat:iEndLat]
    lons = field.grid.get_coords(LON_INDEX, ESMF.StaggerLoc.CORNER)[iBegLon:iEndLon, iBegLat:iEndLat]
    n0, n1 = lats.shape
    numPoints = n0 * n1
    pt.SetNumberOfPoints(numPoints)
    sg.SetDimensions(1, n0, n1)

    # fill in the points
    k = 0
    for i1 in range(n1):
        for i0 in range(n0):
            x = radius*math.cos(lats[i0, i1] * math.pi/180.) * math.cos(lons[i0, i1] * math.pi/180.)
            y = radius*math.cos(lats[i0, i1] * math.pi/180.) * math.sin(lons[i0, i1] * math.pi/180.)
            z = radius*math.sin(lats[i0, i1] * math.pi/180.)
            pt.SetPoint(k, x, y, z)
            k += 1
    # show edges
    ed = vtk.vtkExtractEdges()
    et = vtk.vtkTubeFilter()
    et.SetRadius(0.01)
    em = vtk.vtkPolyDataMapper()
    ea = vtk.vtkActor()
    ed.SetInputData(sg)
    et.SetInputConnection(ed.GetOutputPort())
    em.SetInputConnection(et.GetOutputPort())
    ea.SetMapper(em)
    actors = [ac, ea]
    return actors, mp, sg, pt

def render(actors):
    import vtk
    renderer = vtk.vtkRenderer()
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.AddRenderer(renderer)    

    for ac in actors:
        renderer.AddActor(ac)
    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)
    renderer.SetBackground(.0, .0, .0)
    renderWindow.Render()
    renderWindowInteractor.Start()

# set the src/dst grid dimensions
srcPointDims = (100 + 1, 200 + 1)
dstPointDims = (100 + 1, 200 + 1)

srcCellDims = (srcPointDims[0] - 1, srcPointDims[1] - 1)
dstCellDims = (dstPointDims[0] - 1, dstPointDims[1] - 1)

# 2d latitude and longitude coordinates (rotated grid)
srcLatsPrime = numpy.linspace(-90., 90, srcPointDims[LAT_INDEX])
srcLonsPrime = numpy.linspace(-180., 180, srcPointDims[LON_INDEX])

delta_lat, delta_lon = 10.0, 10.0 #30.0, 20.0
srcLats2D, srcLons2D = createRotatedPoleCoords(srcLatsPrime, srcLonsPrime, delta_lat, delta_lon)

# target grid is lat-lon
dstLatsPrime = numpy.linspace(-89., 89., dstPointDims[LAT_INDEX]) #numpy.linspace(-90., 90, dstPointDims[0])
dstLonsPrime = numpy.linspace(-179., 179., dstPointDims[LON_INDEX]) #numpy.linspace(-180., 180, dstPointDims[1])
dstLats2D, dstLons2D = createRotatedPoleCoords(dstLatsPrime, dstLonsPrime, 0.0, 0.0)

# create the ESMF src/dst grids
srcGrid = createGrid(srcPointDims, srcLats2D, srcLons2D)
dstGrid = createGrid(dstPointDims, dstLats2D, dstLons2D)

# create the ESMF src/dst cell centered fields
srcLats2DCell = 0.25*(srcLats2D[:-1, :-1] + srcLats2D[1:, :-1] + srcLats2D[1:, 1:] + srcLats2D[:-1, 1:])
srcLons2DCell = 0.25*(srcLons2D[:-1, :-1] + srcLons2D[1:, :-1] + srcLons2D[1:, 1:] + srcLons2D[:-1, 1:])
srcData = srcLons2DCell # some arbitrary field
srcField = createField(srcGrid, 'src', srcData)
dstData = numpy.zeros(dstCellDims, numpy.float64) # initializing the dst field to zero
dstField = createField(dstGrid, 'dst', dstData)

# set the field initially to some bad values
dstField.data[...] = -1000.0

plotField(srcField, 'g-')

# compute the interpolation weights
regrid = ESMF.Regrid(srcfield=srcField, dstfield=dstField,
                     regrid_method=ESMF.api.constants.RegridMethod.CONSERVE,
                     unmapped_action=ESMF.api.constants.UnmappedAction.IGNORE)
regrid(srcField, dstField)
print('check sum: {}'.format(dstField.data.sum()))

print(dstField.data)

plotField(dstField, 'r-')

# 3d viz with VTK
#dst_pipeline = plotGrid3d(dstField, color=(0.5, 0., 0.), radius=1.05)
#src_pipeline = plotGrid3d(srcField, color=(0., 0.5, 0.), radius=0.95)
#render(src_pipeline[0] + dst_pipeline[0])

pylab.show()




