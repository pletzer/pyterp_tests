import math
import numpy

def getCellAreas(lats, lons):

    areas1 = numpy.zeros(lats.shape, numpy.float64)
    areas2 = numpy.zeros(lats.shape, numpy.float64)
    
    v1x = lats[1:, :-1] - lats[:-1, :-1]
    v1y = lons[1:, :-1] - lons[:-1, :-1]

    v2x = lats[1:, 1:] - lats[1:, :-1]
    v2y = lons[1:, 1:] - lons[1:, :-1]

    v3x = lats[:-1, 1:] - lats[1:, 1:]
    v3y = lons[:-1, 1:] - lons[1:, 1:]

    v4x = lats[:-1, :-1] - lats[:-1, 1:]
    v4y = lons[:-1, :-1] - lons[:-1, 1:]

    areas1[:-1, :-1] = v1x*v2y - v1y*v2x
    areas2[:-1, :-1] = v3x*v4y - v3y*v4x

    return areas1, areas2



def createCoords(latsPrime, lonsPrime, **kw):
    """
    Create coordinates from axes
    @param latsPrime latitude logical axis
    @param lonsPrime longitude logical axis
    @return curvilinear latitudes, longitudes and data
    """
    delta_lat = kw['delta_lat']
    delta_lon = kw['delta_lon']

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

    # fix the dateline issue. Cells that have negative area must be fixed

    return lats, lons

def createPointData(lats, lons):
    """
    Create data from curvilinear coordinates
    @param lats 2D latitude data
    @param lons 2D longitude data
    @return data
    """
    nj, ni = lats.shape
    data = numpy.zeros(lats.shape, numpy.float64)
    for j in range(nj):
        for i in range(ni):
            # arbitrary function
            data[j, i] = math.sin(2*math.pi*lons[j, i]/180.)*numpy.cos(math.pi*lats[j, i]/180.)

    return data