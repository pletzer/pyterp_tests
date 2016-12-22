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
    rr = numpy.outer(rho, numpy.zeros((nt,), numpy.float64))
    tt = numpy.outer(numpy.zeros((nr,), numpy.float64), the)
    radius = kw.get('radius', 1.0)
    xx = radius * numpy.cos(tt)
    yy = radius * numpy.sin(tt)

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
    xxCells = numpy.zeros((njM1, niM1, 4), numpy.float64)
    data = numpy.zeros((njM1, niM1), numpy.float64)
    for j in range(njM1):
        for i in range(niM1):
            yyCells[j, i, :] = yy[j + 0, i + 0], \
                                yy[j + 0, i + 1], \
                                yy[j + 1, i + 1], \
                                yy[j + 1, i + 0]
            xxCells[j, i, :] = xx[j + 0, i + 0], \
                                xx[j + 0, i + 1], \
                                xx[j + 1, i + 1], \
                                xx[j + 1, i + 0]
            # mid point
            midY = 0.25*yyCells.sum()
            midX = 0.25*xxCells.sum()

            # arbitrary function
            data[j, i] = math.sin(2*math.pi*midX/1.0)*numpy.cos(math.pi*midY/2.0)

    return yyCells, xxCells, data