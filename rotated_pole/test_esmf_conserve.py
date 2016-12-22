import ESMF
import numpy

ESMF.Manager(debug=True)

def create2DLats(nodeDims, latAxis):
    """
    Create a 2d latitude coordinate from a 1d axis
    @param nodeDims number of latitude/longitude nodes
    @param latAxis 1D array of latitude values
    @return 2D lat coordinate
    """
    return numpy.outer(latAxis, numpy.ones((nodeDims[1],)))

def create2DLons(nodeDims, lonAxis):
    """
    Create a 2d longitude coordinate from a 1d axis
    @param nodeDims number of latitude/longitude nodes
    @param lonAxis 1D array of longitude values
    @return 2D lon coordinate
    """
    return numpy.outer(numpy.ones((nodeDims[0],)), lonAxis)

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

# set the src/dst grid dimensions
srcPointDims = (5 + 1, 10 + 1) #(10 + 1, 20 + 1) # (100 + 1, 200 + 1) works!
dstPointDims = (5 + 1, 10 + 1)

srcCellDims = (srcPointDims[0] - 1, srcPointDims[1] - 1)
dstCellDims = (dstPointDims[0] - 1, dstPointDims[1] - 1)

# 2d latitude and longitude coordinates (uniform grids for simplicity)
srcLats2D = create2DLats(srcPointDims, numpy.linspace(-90., 90, srcPointDims[0]))
srcLons2D = create2DLons(srcPointDims, numpy.linspace(-180., 180, srcPointDims[1]))
dstLats2D = create2DLats(dstPointDims, numpy.linspace(-90., 90, dstPointDims[0]))
dstLons2D = create2DLons(dstPointDims, numpy.linspace(-180., 180, dstPointDims[1]))

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

# compute the interpolation weights
regrid = ESMF.Regrid(srcfield=srcField, dstfield=dstField,
                     regrid_method=ESMF.api.constants.RegridMethod.CONSERVE,
                     unmapped_action=ESMF.api.constants.UnmappedAction.ERROR)
#regrid
regrid(srcField, dstField)

print(dstField.data)



