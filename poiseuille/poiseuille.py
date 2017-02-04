#!/usr/bin/env python

"""
Compute the Poiseuille flow inside a pipe of radius 1
as a 2-form (face field)
The grid is Cartesian
"""

import numpy
from math import sqrt, asin, pi
from ctypes import c_int, c_double, c_void_p, CDLL, byref, POINTER
import operator
import nefaceCfg
import os
import sys
import NfError
import NfTriangulator

class Poiseuille:

    def __init__(self, nx, ny, nz, el):
        """
        Constructor
        @param nx number of x cells
        @param ny number of y cells
        @param nx number of z cells
        @param el pipe length
        """

        # open the library
        solib = os.path.normpath(nefaceCfg.libdir + '/' + \
                                     nefaceCfg.libName + \
                                     '.' + nefaceCfg.sharedLibSuffix)
        self.nfc = CDLL(solib)

        # radius of pipe is one
        self.dx = 2.0/float(nx)
        self.dy = 2.0/float(ny)
        self.dz = el/float(nz)
        # nodal coordinates
        self.x = numpy.array([-1.0 + self.dx*i for i in range(nx + 1)])
        self.y = numpy.array([-1.0 + self.dy*j for j in range(ny + 1)])
        self.z = numpy.array([0.0 + self.dz*k for k in range(nz + 1)])

        self.fluxes = numpy.zeros( (nx+1, ny+1, nz+1, 3), numpy.float64 )
        self.flows = numpy.zeros( (nx+1, ny+1, nz+1, 3), numpy.float64 )
        for i in range(nx + 1):
            for j in range(ny + 1):
                fz = self.computeFlowInCell(i, j)
                for k in range(nz + 1):
                    self.fluxes[i, j, k, 0] = 0.0
                    self.fluxes[i, j, k, 1] = 0.0
                    self.fluxes[i, j, k, 2] = fz/(self.dx * self.dy)
                    self.flows[i, j, k, 0] = 0.0
                    self.flows[i, j, k, 1] = 0.0
                    self.flows[i, j, k, 2] = fz

    def saveSurfaceProjection(self, filename, nodes, connectivity, values):
        """
        """
        import tables

        f = tables.openFile(filename, 'w')

        # mesh
        g = f.createGroup('/', name = 'triangleMesh')
        g._v_attrs['vsType'] = 'mesh'
        g._v_attrs['vsKind'] = 'unstructured'
        g._v_attrs['vsPoints'] = '/points'
        g._v_attrs['vsTriangles'] = 'triangles'

        # points
        dPoints = f.createArray('/', 
                                name = 'points', 
                                object = numpy.array(nodes))

        # cells, note connectivity is expected to be int32
        dCells = f.createArray('/triangleMesh', 
                               name = 'triangles', 
                               object = numpy.array(connectivity, numpy.int32))

        # scalar, zonal data of integrated fluxes
        vn = 'flow'
        dVals = f.createArray('/', 
                              name = vn, 
                              object = numpy.array(values, numpy.float64))
        dVals._v_attrs['vsType'] = 'variable'
        dVals._v_attrs['vsMesh'] = '/triangleMesh'
        dVals._v_attrs['vsCentering'] = 'zonal'

        # average field values
        def getArea(ijk):
            d0 = numpy.array(nodes[ijk[0]])
            d1 = numpy.array(nodes[ijk[1]])
            d2 = numpy.array(nodes[ijk[2]])
            d1 -= d0
            d2 -= d0
            crossProd = numpy.cross(d1, d2)
            return 0.5 * numpy.sqrt(numpy.dot(crossProd, crossProd))

        ncells = len(connectivity)
        avgField = [values[icell]/max(1.e-10, getArea(connectivity[icell])) for icell \
                                              in range(ncells)]
        vn = 'areaAveragedField'
        dVals2 = f.createArray('/', 
                              name = vn, 
                              object = numpy.array(avgField, numpy.float64))
        dVals2._v_attrs['vsType'] = 'variable'
        dVals2._v_attrs['vsMesh'] = '/triangleMesh'
        dVals2._v_attrs['vsCentering'] = 'zonal'

        # done
        f.close()          

    def save(self, filename):
        """
        Save face data to VizSchema compliant HDF5 file
        """
        import tables
        
        h5file = tables.openFile(filename, mode = "w")

        ggg = h5file.createGroup("/", "globalGridGlobal")
        ggg._f_setAttr("vsType", "mesh")
        ggg._f_setAttr("vsKind", "uniform")
        nx, ny, nz = len(self.x) - 1, len(self.y) - 1, len(self.z) - 1
        ggg._f_setAttr("vsNumCells", numpy.array((nx, ny, nz), 
                                                 numpy.int32))
        ggg._f_setAttr("vsLowerBounds", numpy.array((self.x[0], 
                                                     self.y[0], 
                                                     self.z[0])))
        ggg._f_setAttr("vsUpperBounds",  numpy.array((self.x[-1], 
                                                      self.y[-1], 
                                                      self.z[-1])))
        ggg._f_setAttr("vsStartCell", numpy.array((0, 0, 0)))

        gggl = h5file.createGroup("/", "globalGridGlobalLimits")
        gggl._f_setAttr("vsKind", "Cartesian")
        gggl._f_setAttr("vsType", "limits")
        gggl._f_setAttr("vsLowerBounds", numpy.array((self.x[0], self.y[0], self.z[0])))
        gggl._f_setAttr("vsUpperBounds", numpy.array((self.x[-1], 
                                                      self.y[-1], 
                                                      self.z[-1])))
        
        d = h5file.createArray("/", "fluxes", self.fluxes)
        d.attrs.vsType = "variable"
        d.attrs.vsMesh = "globalGridGlobal"
        d.attrs.vsLimits = "globalGridGlobalLimits"
        d.attrs.vsCentering = "face"

        d2 = h5file.createArray("/", "flows", self.flows)
        d2.attrs.vsType = "variable"
        d2.attrs.vsMesh = "globalGridGlobal"
        d2.attrs.vsLimits = "globalGridGlobalLimits"
        d2.attrs.vsCentering = "face"

        dv = h5file.createGroup("/", "derivedVariables")
        dv._f_setAttr("vsType", "vsVars")
        dv._f_setAttr("Flux", "{fluxes_0, fluxes_1, fluxes_2}")
        dv._f_setAttr("Flow", "{flows_0, flows_1, flows_2}")
        
        h5file.close()
        

    def projectOntoSurface(self, expressions, nu, nv):
        """
        Project field onto a surface
        @param expressions [x(u,v), y(u, v), z(u, v)]
        @param number of u cells
        @param number of v cells
        @return nodes, connectivity, values
        """
        xmins = numpy.array( (self.x[0], self.y[0], self.z[0]) )
        xmaxs = numpy.array( (self.x[-1], self.y[-1], self.z[-1]) )
        ns = numpy.array((len(self.x), len(self.y), len(self.z)), 
                         numpy.int32)
        tr = NfTriangulator.Triangulator(expressions,
                                         xmins = xmins, 
                                         xmaxs = xmaxs, 
                                         ns = ns, 
                                         nUCells = nu,
                                         nVCells = nv)
        tr.triangulate()
        triangles = tr.getTriangles()
        nodes = tr.getPoints()
        connectivity = tr.getConnectivity()

        n = len(triangles)
        values = numpy.zeros( (n,), numpy.float64 )

        npdims = 3
        ndims = 3

        # handles
        mesh_t = c_void_p(0)
        field_t = c_void_p(0)
        nodeDims = [n for n in ns]
        cellDims_t = (c_int * npdims)()
        cellDims_t[:] = [n-1 for n in ns]
        npoints = reduce(operator.mul, nodeDims)

        # create uniform mesh
        ier = self.nfc.NfStructuredGridNew(byref(mesh_t), 
                                           npdims, ndims)
        NfError.fatal(ier, self.nfc.NfStructuredGridError, __file__, 
                      sys._getframe().f_lineno)
        ier = self.nfc.NfStructuredGridSetCellDimensions(byref(mesh_t), 
                                                         cellDims_t)
        NfError.fatal(ier, self.nfc.NfStructuredGridError, __file__, 
                      sys._getframe().f_lineno)
        points = numpy.zeros( (npoints, ndims), numpy.float64 )
        inode = 0
        for i in range(nodeDims[0]):
            for j in range(nodeDims[1]):
                for k in range(nodeDims[2]):
                    points[inode, :] = self.x[i], self.y[j], self.z[k]
                    inode += 1
                    
        points_t = points.ctypes.data_as(POINTER(c_double))
        ier = self.nfc.NfStructuredGridSetVertices(byref(mesh_t), points_t)
        NfError.fatal(ier, self.nfc.NfStructuredGridError, __file__, 
                      sys._getframe().f_lineno)

        # create uniform field
        ier = self.nfc.NfStructuredFieldNew(byref(field_t), 
                                      "face_field",
                                      mesh_t)
        NfError.fatal(ier, self.nfc.NfStructuredFieldError, __file__, 
                      sys._getframe().f_lineno)
        
        xFace_t = (c_int*npdims)(0, -1, -1)
        yFace_t = (c_int*npdims)(-1, 0, -1)
        zFace_t = (c_int*npdims)(-1, -1, 0)
        # fluxes across faces, make sure we're getting some 
        # fresh contiguous arrays
        xData = self.flows[..., 0].copy()
        yData = self.flows[..., 1].copy()
        zData = self.flows[..., 2].copy()
        xData_t = xData.ctypes.data_as(POINTER(c_double))
        yData_t = yData.ctypes.data_as(POINTER(c_double))
        zData_t = zData.ctypes.data_as(POINTER(c_double))
        ier = self.nfc.NfStructuredFieldSetPointer(byref(field_t), 
                                                   xFace_t, xData_t)
        NfError.fatal(ier, self.nfc.NfStructuredFieldError, __file__, 
                      sys._getframe().f_lineno)
        ier = self.nfc.NfStructuredFieldSetPointer(byref(field_t), 
                                                   yFace_t, yData_t)
        NfError.fatal(ier, self.nfc.NfStructuredFieldError, __file__, 
                      sys._getframe().f_lineno)
        ier = self.nfc.NfStructuredFieldSetPointer(byref(field_t), 
                                                   zFace_t, zData_t)
        NfError.fatal(ier, self.nfc.NfStructuredFieldError, __file__, 
                      sys._getframe().f_lineno)
        #self.nfc.NfStructuredFieldPrint(byref(field_t)) #####
        
        # iterate over triangles
        val_t = c_double(0)
        p0_t = (c_double * ndims)()
        p1_t = (c_double * ndims)()
        p2_t = (c_double * ndims)()
        i0_t = (c_double * ndims)()
        i1_t = (c_double * ndims)()
        i2_t = (c_double * ndims)()

        for itri in range(len(triangles)):

            p0_t[:] = triangles[itri][0]
            p1_t[:] = triangles[itri][1]
            p2_t[:] = triangles[itri][2]

            #
            # find the vertex coordinates in index space
            #

            i0_t[:] = [d/2 for d in nodeDims]
            ier = self.nfc.NfStructuredGridFindCell(byref(mesh_t), 
                                                    p0_t, i0_t)
            if ier != 0:
                print 'Failed to find index of point p0 for triangle %d p0 = %f %f %f' % \
                    (itri, p0_t[0], p0_t[1], p0_t[2])
                print 'Point outside of domain?'
                NfError.bad(ier, self.nfc.NfStructuredGridError, __file__, 
                            sys._getframe().f_lineno)

            i1_t[:] = [d/2 for d in nodeDims]
            ier = self.nfc.NfStructuredGridFindCell(byref(mesh_t), 
                                                    p1_t, i1_t)
            if ier != 0:
                print 'Failed to find index of point p1 for triangle %d p1 = %f %f %f' % \
                    (itri, p1_t[0], p1_t[1], p1_t[2])
                print 'Point outside of domain?'
                NfError.bad(ier, self.nfc.NfStructuredGridError, __file__, 
                            sys._getframe().f_lineno)

            i2_t[:] = [d/2 for d in nodeDims]
            ier = self.nfc.NfStructuredGridFindCell(byref(mesh_t), 
                                                    p2_t, i2_t)
            if ier != 0:
                print 'Failed to find index of point p2 for triangle %d p2 = %f %f %f' % \
                    (itri, p2_t[0], p2_t[1], p2_t[2])
                print 'Point outside of domain?'
                NfError.bad(ier, self.nfc.NfStructuredGridError, __file__, 
                            sys._getframe().f_lineno)

            # now project
            ier = self.nfc.NfStructuredFieldProjectToTriangle(byref(field_t), 
                                                              i0_t, 
                                                              i1_t, 
                                                              i2_t,
                                                              byref(val_t))
            NfError.fatal(ier, self.nfc.NfStructuredFieldError, __file__, 
                          sys._getframe().f_lineno)

            values[itri] = val_t.value

        # clean up
        ier = self.nfc.NfStructuredFieldDelete(byref(field_t))
        NfError.fatal(ier, self.nfc.NfStructuredFieldError, __file__, 
                      sys._getframe().f_lineno)
        ier = self.nfc.NfStructuredGridDelete(byref(mesh_t))
        NfError.fatal(ier, self.nfc.NfStructuredGridError, __file__, 
                      sys._getframe().f_lineno)

        return nodes, connectivity, values
        
        


    def computeTotalFlow(self):
        """
        Compute total flow in pipe
        """
        return numpy.sum(self.fluxes[:,:,0,2].flat) * self.dx * self.dy

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
def test0():
    p = Poiseuille(3, 3, 2, 2.0)
    print p.computeFlowInCutCircleCell(x0 = -1./3., x1 = 1./3., y0 = sqrt(1. - 1./9.))

def test1():
    el = 4.0
    p = Poiseuille(5, 6, 3, el)
    totFlow = p.computeTotalFlow()
    exactTotalFlow = 2*pi*0.25
    print 'totFlow = ', totFlow, ' exact = ', exactTotalFlow, ' error = ', totFlow - exactTotalFlow
    p.save('poiseuille.vsh5')
    nodes, connectivity, values = p.projectOntoSurface( \
        ['0.1+0.8*u*cos(2*pi*v)', '0.2+0.8*u*sin(2*pi*v)', 'zmin+zlen*v'], \
        10, 20)
    print 'projected flow = ', numpy.sum(values)
    p.saveSurfaceProjection('poiseuilleSurf.vsh5', nodes, connectivity, values)

if __name__ == '__main__': 
    test1()
