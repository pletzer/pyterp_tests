import math
import numpy

def createCoordAndData(latsPrime, lonsPrime, **kw):
    """
    Create coordinates and data from axes
    @param latsPrime latitude logical axis
    @param lonsPrime longitude logical axis
    @return curvilinear latitudes, longitudes and data
    """
    delta_theta = kw.get('delta_theta', 0.0)
    delta_lambda = kw.get('delta_lambda', 0.0)

    nj, ni = len(latsPrime), len(lonsPrime)
    lats = numpy.zeros((nj, ni,), numpy.float64)
    lons = numpy.zeros((nj, ni,), numpy.float64)
    data = numpy.zeros((nj, ni,), numpy.float64)
    for j in range(nj):
        the = math.pi * (latsPrime[j] - delta_theta) / 180.
        cos_the = math.cos(the)
        sin_the = math.sin(the)
        for i in range(ni):
            lam = math.pi * (lonsPrime[i] - delta_lambda) / 180.
            cos_lam = math.cos(lam)
            sin_lam = math.sin(lam)
            x = cos_the * cos_lam
            y = cos_the * sin_lam
            z = sin_the
            rho = math.sqrt(x*x + y*y)
            lats[j, i] = 180. * math.atan2(z, rho) / math.pi
            lons[j, i] = 180. * math.atan2(y, x) / math.pi
            # arbitrary function
            data[j, i] = math.sin(2*math.pi*lons[j, i]/180.)*numpy.cos(math.pi*lats[j, i]/180.)

    return lats, lons, data
