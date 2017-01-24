import ESMF
import numpy
import math
from matplotlib import pylab

ESMF.Manager(debug=True)

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
    lats = numpy.zeros((nj, ni,), numpy.float64)
    lons = numpy.zeros((nj, ni,), numpy.float64)

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

            lats[j, i] = 180. * math.asin(xyz[2]) / math.pi
            lons[j, i] = 180. * math.atan2(xyz[1], xyz[0]) / math.pi

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
    latIndex, lonIndex = 0, 1

    grid.add_coords(staggerloc=ESMF.StaggerLoc.CORNER, coord_dim=latIndex)
    coordLatPoint = grid.get_coords(coord_dim=latIndex, staggerloc=ESMF.StaggerLoc.CORNER)

    grid.add_coords(staggerloc=ESMF.StaggerLoc.CORNER, coord_dim=lonIndex)
    coordLonPoint = grid.get_coords(coord_dim=lonIndex, staggerloc=ESMF.StaggerLoc.CORNER)

    iBegLat = grid.lower_bounds[ESMF.StaggerLoc.CORNER][latIndex]
    iEndLat = grid.upper_bounds[ESMF.StaggerLoc.CORNER][latIndex]
    iBegLon = grid.lower_bounds[ESMF.StaggerLoc.CORNER][lonIndex]
    iEndLon = grid.upper_bounds[ESMF.StaggerLoc.CORNER][lonIndex]

    coordLatPoint[...] = lats2D[iBegLat:iEndLat, iBegLon:iEndLon]
    coordLonPoint[...] = lons2D[iBegLat:iEndLat, iBegLon:iEndLon]

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
    latIndex, lonIndex = 0, 1
    iBegLat = field.grid.lower_bounds[ESMF.StaggerLoc.CORNER][latIndex]
    iEndLat = field.grid.upper_bounds[ESMF.StaggerLoc.CORNER][latIndex]
    iBegLon = field.grid.lower_bounds[ESMF.StaggerLoc.CORNER][lonIndex]
    iEndLon = field.grid.upper_bounds[ESMF.StaggerLoc.CORNER][lonIndex]

    lats = field.grid.get_coords(latIndex, ESMF.StaggerLoc.CORNER)[iBegLat:iEndLat, iBegLon:iEndLon]
    lons = field.grid.get_coords(lonIndex, ESMF.StaggerLoc.CORNER)[iBegLat:iEndLat, iBegLon:iEndLon]

    # plot the grid
    nc = 20 # number of curves
    nj, ni = lats.shape
    stepj, stepi = max(1, nj//nc), max(1, ni//nc)
    for j in range(0, nj, stepj):
        pylab.plot(lons[j, :], lats[j, :], linetype)
    for i in range(0, ni, stepi):
        pylab.plot(lons[:, i], lats[:, i], linetype)

    pylab.pcolormesh(lons, lats, field.data[0,...])
    #pylab.show()

# set the src/dst grid dimensions
srcPointDims = (100 + 1, 200 + 1)
dstPointDims = (5 + 1, 10 + 1)

srcCellDims = (srcPointDims[0] - 1, srcPointDims[1] - 1)
dstCellDims = (dstPointDims[0] - 1, dstPointDims[1] - 1)

# 2d latitude and longitude coordinates (rotated grid)
srcLatsPrime = numpy.linspace(-90., 90, srcPointDims[0])
srcLonsPrime = numpy.linspace(-180., 180, srcPointDims[1])

delta_lat, delta_lon = 10.0, 10.0 #30.0, 20.0
srcLats2D, srcLons2D = createRotatedPoleCoords(srcLatsPrime, srcLonsPrime, delta_lat, delta_lon)

# target grid is lat-lon
dstLatsPrime = numpy.linspace(-20., 20, dstPointDims[0]) #numpy.linspace(-90., 90, dstPointDims[0])
dstLonsPrime = numpy.linspace(-110., 110, dstPointDims[1]) #numpy.linspace(-180., 180, dstPointDims[1])
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
#regrid
regrid(srcField, dstField)

print(dstField.data)

plotField(dstField, 'r-')
pylab.show()




