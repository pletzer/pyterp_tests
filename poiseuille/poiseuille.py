#!/usr/bin/env python

"""
Compute the Poiseuille flow inside a quater of a pipe of radius 1
as a 2-form (face field)
The grid is Cartesian
"""

import numpy
from math import sqrt, asin, pi
from ctypes import c_int, c_double, c_void_p, CDLL, byref, POINTER
import operator
import os
import sys
import argparse

class Poiseuille:

    def __init__(self, nx1, ny1):
        """
        Constructor
        @param nx number of x points
        @param ny number of y points
        """
        # radius of pipe is one
        self.dx = 1.0/float(nx1 - 1)
        self.dy = 1.0/float(ny1 - 1)
        areaElement = self.dx * self.dy

        # nodal coordinates, domain is [0, 1] x [0, 1]
        self.x = numpy.array([0.0 + self.dx*i for i in range(nx1)])
        self.y = numpy.array([0.0 + self.dy*j for j in range(ny1)])
        self.yy = numpy.zeros( (nx1, ny1,), numpy.float64 )
        self.xx = numpy.zeros( (nx1, ny1,), numpy.float64 )
        self.field = numpy.zeros( (nx1, ny1,), numpy.float64 )
        self.flux = numpy.zeros( (nx1 - 1, ny1 - 1,), numpy.float64 )

        for i in range(nx1):
            for j in range(ny1):
                self.yy[i, j] = self.y[j]
                self.xx[i, j] = self.x[i]
                self.field[i, j] = 1.0 - self.x[i]**2 - self.y[j]**2

        for i in range(nx1 - 1):
            for j in range(ny1 - 1):
                # most conservative interpolation schemes want 
                # fluxes
                self.flux[i, j] = self.computeFlowInCell(i, j) / areaElement

    def save(self, filename):
        """
        Save flow in netcdf file
        """
        import iris

        yBounds, yMid = self.computeBounds(self.yy)
        xBounds, xMid = self.computeBounds(self.xx)
        yCoordMid = iris.coords.AuxCoord(yMid, var_name='yMid', standard_name='latitude', units='1', bounds=xBounds)
        xCoordMid = iris.coords.AuxCoord(xMid, var_name='xMid', standard_name='longitude', units='1', bounds=yBounds)
        yCoord = iris.coords.AuxCoord(self.yy, var_name='yy', standard_name='latitude', units='1')
        xCoord = iris.coords.AuxCoord(self.xx, var_name='xx', standard_name='longitude', units='1')

        cubeCell = iris.cube.Cube(self.flux, var_name='flux', standard_name='stratiform_precipitation_flux', units='kg/m^2')
        cubeCell.add_aux_coord(yCoordMid, data_dims=(0, 1))
        cubeCell.add_aux_coord(xCoordMid, data_dims=(0, 1))

        cubePoint = iris.cube.Cube(self.field, var_name='pointData', standard_name='stratiform_precipitation_flux', units='kg/m^2')
        cubePoint.add_aux_coord(yCoord, data_dims=(0, 1))
        cubePoint.add_aux_coord(xCoord, data_dims=(0, 1))

        iris.save([cubeCell, cubePoint], filename)

    def computeBounds(self, array2d):
        """
        From the point coordinates, compute the bounds and mid point values
        @param array2d coordinate as a 2d array
        @return bounds, mid points
        """
        # number of cells
        nx, ny = array2d.shape[0] - 1, array2d.shape[1] - 1
        coordBounds = numpy.zeros((nx, ny, 4), numpy.float64)
        coordMid = numpy.zeros((nx, ny), numpy.float64)
        for i in range(nx):
            for j in range(ny):
                coordBounds[i, j, 0] = array2d[i  , j  ]
                coordBounds[i, j, 1] = array2d[i+1, j  ]
                coordBounds[i, j, 2] = array2d[i+1, j+1]
                coordBounds[i, j, 3] = array2d[i  , j+1]
                coordMid[i, j] = 0.25*coordBounds[i, j, :].sum()
        return coordBounds, coordMid

    def computeTotalFlow(self):
        """
        Compute total flow in pipe
        """
        return numpy.sum(self.flux[:,:].flat) * self.dx * self.dy

    def computeFlowInCell(self, i, j):
        """
        Compute the flow in a cell
        @param i low x-index
        @param j low y-index
        """
        # cell boundary
        x0 = self.x[i]
        x1 = x0 + self.dx
        y0 = self.y[j]
        y1 = y0 + self.dy        
        numOutsideNodes = self.getNumberOfOutsideNodes(i, j)
        if numOutsideNodes == 0:
            return self.computeFlowInFullCell(x0, x1, y0, y1)
        elif numOutsideNodes == 1:
            return self.computeFlowInPartialCell1(x0, x1, y0, y1)
        elif numOutsideNodes == 2:
            return self.computeFlowInPartialCell2(x0, x1, y0, y1)
        elif numOutsideNodes == 3:
            return self.computeFlowInPartialCell3(x0, x1, y0, y1)
        # invalid cell
        return 0.0

    def computeFlowInFullCell(self, x0, x1, y0, y1):
        """
        Compute integrated flux in a full cell
        @param x0 lower x-bound
        @param x1 upper x-bound
        @param y0 lower y-bound
        @param y1 upper y-bound
        @return flow
        """
        dx = x1 - x0
        dy = y1 - y0
        res = dx*dy
        res -= (x1**3 - x0**3)*dy/3.0
        res -= dx*(y1**3 - y0**3)/3.0
        return res


    def computeFlowInCutCircleCell(self, x0, x1, y0):
        """
        Compute the flow through a disk cut along x and y
        @param x0 lower integration endpoint
        @param x1 upper integration endpoint
        @param y0 lower y-integration endpoint (upper is the circle)
        """
        res = (-(x0*sqrt(1 - x0**2)) - \
                    sqrt(1 - x0**2)*(x0/4. - x0**3/2.) + 3*x0*y0 - \
                    x0**3*y0 - x0*y0**3 - (3*asin(x0))/4.)/3. + \
                    (x1*sqrt(1 - x1**2) + \
                         sqrt(1 - x1**2)*(x1/4. - x1**3/2.) - 3*x1*y0 + \
                         x1**3*y0 + x1*y0**3 + (3*asin(x1))/4.)/3.
        return res
    

    def computeFlowInPartialCell1(self, x0, x1, y0, y1):

        # 1 node is invalid

        res = 0.0

        if x0**2 + y0**2 > 1.0:
            x2 = self.findXEdgeIntersect(x0, y0)
            y2 = self.findYEdgeIntersect(y0, x0)
            res += self.computeFlowInCutCircleCell(-x2, -x0, -y2)
            res += self.computeFlowInFullCell(x2, x1, y0, y1)
            res += self.computeFlowInFullCell(x0, x2, y2, y1)
        elif x1**2 + y0**2 > 1.0:
            x2 = self.findXEdgeIntersect(x1, y0)
            y2 = self.findYEdgeIntersect(y0, x1)
            res += self.computeFlowInCutCircleCell(x2, x1, y2)
            res += self.computeFlowInFullCell(x0, x2, y0, y1)
            res += self.computeFlowInFullCell(x2, x1, y2, y1)
        elif x1**2 + y1**2 > 1.0:
            x2 = self.findXEdgeIntersect(x1, y1)
            y2 = self.findYEdgeIntersect(y1, x1)
            res += self.computeFlowInCutCircleCell(x2, x1, y2)
            res += self.computeFlowInFullCell(x0, x2, y0, y1)
            res += self.computeFlowInFullCell(x2, x1, y0, y2)
        elif x0**2 + y1**2 > 1.0:
            x2 = self.findXEdgeIntersect(x0, y1)
            y2 = self.findYEdgeIntersect(y1, x0)
            res -= self.computeFlowInCutCircleCell(x2, x0, y2)
            res += self.computeFlowInFullCell(x2, x1, y0, y1)
            res += self.computeFlowInFullCell(x0, x2, y0, y2)
       
        return res

    def computeFlowInPartialCell2(self, x0, x1, y0, y1):

        # 2 nodes are invalid
        
        res = 0.0
        if x0**2 + max(y0**2, y1**2) < 1.0:
            # left is valid
            x2 = self.findXEdgeIntersect(x0, y0)
            x3 = self.findXEdgeIntersect(x0, y1)
            res += self.computeFlowInFullCell(x0, min(x2, x3), y0, y1)
            res += self.computeFlowInCutCircleCell(y0, y1, min(x2, x3))
        elif max(x0**2, x1**2) + y0**2 < 1.0:
            # down side is valid
            y2 = self.findYEdgeIntersect(y0, x0)
            y3 = self.findYEdgeIntersect(y0, x1)
            res += self.computeFlowInFullCell(x0, x1, y0, min(y2, y3))
            res += self.computeFlowInCutCircleCell(x0, x1, min(y2, y3))
        elif x1**2 + max(y0**2, y1**2) < 1.0:
            # right side is valid
            x2 = self.findXEdgeIntersect(x0, y0)
            x3 = self.findXEdgeIntersect(x0, y1)
            res += self.computeFlowInFullCell(max(x2, x3), x1, y0, y1)
            res += self.computeFlowInCutCircleCell(y0, y1, min(-x2, -x3))
        elif max(x0**2, x1**2) + y1**2 < 1.0:
            # upper side is valid
            y2 = self.findYEdgeIntersect(y0, x0)
            y3 = self.findYEdgeIntersect(y0, x1)
            res += self.computeFlowInFullCell(x0, x1, max(y2, y3), y1)
            res += self.computeFlowInCutCircleCell(x0, x1, min(-y2, -y3))
            
        return res

    def computeFlowInPartialCell3(self, x0, x1, y0, y1):

        # 3 nodes are invalid

        res = 0.0
        if x0**2 + y0**2 < 1.0:
            # lower left is valid
            x2 = self.findXEdgeIntersect(x0, y0)
            y2 = self.findYEdgeIntersect(y0, x0)
            res += self.computeFlowInCutCircleCell(x0, x2, y0)
        elif x1**2 + y0**2 < 1.0:
            # lower right is valid
            x2 = self.findXEdgeIntersect(x1, y0)
            y2 = self.findYEdgeIntersect(y0, x1)
            res += self.computeFlowInCutCircleCell(x2, x1, y0)
        elif x1**2 + y1**2 < 1.0:
            # upper right is valid
            x2 = self.findXEdgeIntersect(x1, y1)
            y2 = self.findYEdgeIntersect(y1, x1)
            res += self.computeFlowInCutCircleCell(-x1, -x2, -y1)
        elif x0**2 + y1**2 < 1.0:
            # upper left is valid
            x2 = self.findXEdgeIntersect(x0, y1)
            y2 = self.findYEdgeIntersect(y1, x0)
            res += self.computeFlowInCutCircleCell(x0, x2, -y1)

        return res

    def findXEdgeIntersect(self, closeX, y):
        res = sqrt(1.0 - y**2)
        if abs(closeX - res) < abs(closeX + res):
            return res
        return -res
        
    def findYEdgeIntersect(self, closeY, x):
        res = sqrt(1.0 - x**2)
        if abs(closeY - res) < abs(closeY + res):
            return res
        return -res
        
    def getNumberOfOutsideNodes(self, i, j):
        x0 = self.x[i]
        x1 = x0 + self.dx
        y0 = self.y[j]
        y1 = y0 + self.dy
        res = 0
        if x0**2 + y0**2 > 1.0: res +=1
        if x1**2 + y0**2 > 1.0: res +=1
        if x0**2 + y1**2 > 1.0: res +=1
        if x1**2 + y1**2 > 1.0: res +=1
        return res

################################################################################
def main():
    parser = argparse.ArgumentParser(description='Interpolate using ESMF')
    parser.add_argument('--src_file', type=str, dest='src_file', default='src.nc',
                help='Source data file name')
    parser.add_argument('--src_nj', type=int, dest='src_nj', default=11, 
                help='Number of source cells in the y direction')
    parser.add_argument('--src_ni', type=int, dest='src_ni', default=11, 
                help='Number of source cells in the x direction')
    parser.add_argument('--dst_file', type=str, dest='dst_file', default='dst.nc',
                    help='Destination data file name')
    parser.add_argument('--dst_nj', type=int, dest='dst_nj', default=21, 
                help='Number of destination cells in the y direction')
    parser.add_argument('--dst_ni', type=int, dest='dst_ni', default=21, 
                help='Number of destination cells in the x direction')

    args = parser.parse_args()


    # quarter of disc
    exactTotalFlow = 2*pi*0.25 / 4. 

    # src grid
    src = Poiseuille(args.src_nj, args.src_ni)
    totFlow = src.computeTotalFlow()
    print 'src totFlow = ', totFlow, ' exact = ', exactTotalFlow, ' error = ', totFlow - exactTotalFlow
    src.save('src.nc')

    # dst grid
    dst = Poiseuille(args.dst_nj, args.dst_ni)
    totFlow = src.computeTotalFlow()
    print 'src totFlow = ', totFlow, ' exact = ', exactTotalFlow, ' error = ', totFlow - exactTotalFlow
    src.save('dst.nc')

if __name__ == '__main__': 
    main()
