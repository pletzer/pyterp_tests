import math
import numpy

def getCellAreas(yy, xx):

    areas1 = numpy.zeros(yy.shape, numpy.float64)
    areas2 = numpy.zeros(yy.shape, numpy.float64)
    
    v1x = yy[1:, :-1] - yy[:-1, :-1]
    v1y = xx[1:, :-1] - xx[:-1, :-1]

    v2x = yy[1:, 1:] - yy[1:, :-1]
    v2y = xx[1:, 1:] - xx[1:, :-1]

    v3x = yy[:-1, 1:] - yy[1:, 1:]
    v3y = xx[:-1, 1:] - xx[1:, 1:]

    v4x = yy[:-1, :-1] - yy[:-1, 1:]
    v4y = xx[:-1, :-1] - xx[:-1, 1:]

    areas1[:-1, :-1] = v1x*v2y - v1y*v2x
    areas2[:-1, :-1] = v3x*v4y - v3y*v4x

    return areas1, areas2



def createCoords(rho, the, **kw):
    """
    Create coordinates from axes
    @param rho radial axis
    @param the poloidal axis
    @return curvilinear y and x
    """
    nr, nt = len(rho), len(the)
    rr = numpy.outer(rho, numpy.ones((nt,), numpy.float64))
    tt = numpy.outer(numpy.ones((nr,), numpy.float64), the)
    radius = kw.get('radius', 1.0)
    xx = rr * numpy.cos(tt)
    yy = rr * numpy.sin(tt)

    return yy, xx

def createPointData(yy, xx):
    """
    Create nodal data from curvilinear coordinates
    @param yy 2D y data
    @param xx 2D x data
    @return data
    """
    nj, ni = yy.shape
    # arbitrary function
    data = numpy.sin(2*numpy.pi*xx/1.0) * numpy.cos(2*numpy.pi*yy/2.0)
    return data

def createCellData(yy, xx):
    """
    Create zonal data from curvilinear coordinates
    @param yy 2D y data
    @param xx 2D x data
    @return yyCells, xxCells, data
    """
    nj, ni = yy.shape
    njM1 , niM1 = nj - 1, ni - 1

    yyCells = numpy.zeros((njM1, niM1, 4), numpy.float64)
    yyCells[..., 0] = yy[:-1, :-1]
    yyCells[..., 1] = yy[:-1, 1:]
    yyCells[..., 2] = yy[1:, 1:]
    yyCells[..., 3] = yy[1:, :-1]

    xxCells = numpy.zeros((njM1, niM1, 4), numpy.float64)
    xxCells[..., 0] = xx[:-1, :-1]
    xxCells[..., 1] = xx[:-1, 1:]
    xxCells[..., 2] = xx[1:, 1:]
    xxCells[..., 3] = xx[1:, :-1]

    # mid point
    yy = 0.25*yyCells.sum(axis=2)
    xx = 0.25*xxCells.sum(axis=2)

    # arbitrary function
    data = numpy.sin(2*math.pi*xx/1.0)*numpy.cos(math.pi*yy/2.0)

    return yyCells, xxCells, data