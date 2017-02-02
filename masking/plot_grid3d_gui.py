#!/usr/bin/env python

"""
Plot grid in 3D

Requires: NumPy, Iris, VTK and PyQt5

"""
from __future__ import print_function
import math

import iris
import numpy as np
from PyQt5 import QtWidgets, QtCore
import vtk
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor


iris.FUTURE.netcdf_promote = True


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self, srcFile, dstFile):
        super(MainWindow, self).__init__()
        
        # resize window
        self.resize(QtCore.QSize(1000, 600))
        
        # menus
        self.viewMenu = self.menuBar().addMenu("&View")
        
        # settings
        self.src_radius = 0.99
        self.dst_radius = 1.05
        self.resolution = 16
        # TODO: these should depend on number of grid points
        self.src_field_radius = 0.02
        self.dst_field_radius = 0.01
        
        # VTK window
        self._vtkWidget = QVTKRenderWindowInteractor()
        self._vtkWidget._Iren.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())
        self._vtkWidget.Initialize()
        self._vtkWidget.Start()
        self._vtkRenWin = self._vtkWidget.GetRenderWindow()
        self._vtkRen = vtk.vtkRenderer()
        self._vtkRen.SetBackground(0,0,0)
        self._vtkRenWin.AddRenderer(self._vtkRen)
        self.setCentralWidget(self._vtkWidget)
        
        # read data
        self._src_cube = self.readCube(srcFile)
        self._dst_cube = self.readCube(dstFile)
        
        # make vtk points/polydata
        self.field_min = None
        self.field_max = None
        self._src_point_data = self.makePoints(self._src_cube, self.src_radius)
        self._dst_point_data = self.makePoints(self._dst_cube, self.dst_radius)
        
        # create actors for grid
        src_pipeline = self.createGrid(self._src_cube, self._src_point_data[0], color=(0.0, 1.0, 0.0),
                                       radius=self.src_radius, show_mesh_as_surface=True)
        self._src_grid_actors = src_pipeline[0]
        dst_pipeline = self.createGrid(self._dst_cube, self._dst_point_data[0], color=(1.0, 0.0, 0.0),
                                       radius=self.dst_radius)
        self._dst_grid_actors = dst_pipeline[0]
        
        # create actors for field
        self._src_field_actors = self.createField(self._src_point_data[1], self._src_point_data[2],
                                                  self.src_field_radius, scalar_bar=True)
        self._dst_field_actors = self.createField(self._dst_point_data[1], self._dst_point_data[2],
                                                  self.dst_field_radius)
        
        # settings
        dw = self.init_dock_widget()
        dock = QtWidgets.QDockWidget("Display options", self)
        dock.setAllowedAreas(QtCore.Qt.LeftDockWidgetArea | QtCore.Qt.RightDockWidgetArea)
        dock.setWidget(dw)
        self.addDockWidget(QtCore.Qt.RightDockWidgetArea, dock)
        self.viewMenu.addAction(dock.toggleViewAction())
        
        # show the window
        self.setFocus()
    
    def toggle_src_grid(self, state):
        """Toggle src grid visibility."""
        if state == QtCore.Qt.Unchecked:
            for actor in self._src_grid_actors:
                if actor is not None:
                    self._vtkRen.RemoveActor(actor)
        
        else:
            for actor in self._src_grid_actors:
                if actor is not None:
                    self._vtkRen.AddActor(actor)
        
        self._vtkWidget.ReInitialize()
    
    def toggle_dst_grid(self, state):
        """Toggle dst grid visibility."""
        if state == QtCore.Qt.Unchecked:
            for actor in self._dst_grid_actors:
                if actor is not None:
                    self._vtkRen.RemoveActor(actor)
        
        else:
            for actor in self._dst_grid_actors:
                if actor is not None:
                    self._vtkRen.AddActor(actor)
        
        self._vtkWidget.ReInitialize()
    
    def toggle_src_field(self, state):
        """Toggle src field visibility."""
        if state == QtCore.Qt.Unchecked:
            for actor in self._src_field_actors:
                if actor is not None:
                    self._vtkRen.RemoveActor(actor)
        
        else:
            for actor in self._src_field_actors:
                if actor is not None:
                    self._vtkRen.AddActor(actor)
        
        self._vtkWidget.ReInitialize()
    
    def toggle_dst_field(self, state):
        """Toggle dst field visibility."""
        if state == QtCore.Qt.Unchecked:
            for actor in self._dst_field_actors:
                if actor is not None:
                    self._vtkRen.RemoveActor(actor)
        
        else:
            for actor in self._dst_field_actors:
                if actor is not None:
                    self._vtkRen.AddActor(actor)
        
        self._vtkWidget.ReInitialize()
    
    def init_dock_widget(self):
        """Initialise dock widget."""
        # container
        show_widget = QtWidgets.QWidget()
        show_layout = QtWidgets.QVBoxLayout(show_widget)
        
        # src options
        src_group = QtWidgets.QGroupBox("src options", parent=show_widget)
        vbox = QtWidgets.QVBoxLayout()
        # src grid
        src_grid_check = QtWidgets.QCheckBox("Show grid")
        src_grid_check.stateChanged.connect(self.toggle_src_grid)
        src_grid_check.setCheckState(QtCore.Qt.Checked)
        vbox.addWidget(src_grid_check)
        # src field
        src_field_check = QtWidgets.QCheckBox("Show field")
        src_field_check.stateChanged.connect(self.toggle_src_field)
        src_field_check.setCheckState(QtCore.Qt.Checked)
        vbox.addWidget(src_field_check)
        # set group layout
        src_group.setLayout(vbox)
        
        # dst options
        dst_group = QtWidgets.QGroupBox("dst options", parent=show_widget)
        vbox = QtWidgets.QVBoxLayout()
        # dst grid
        dst_grid_check = QtWidgets.QCheckBox("Show grid")
        dst_grid_check.stateChanged.connect(self.toggle_dst_grid)
        dst_grid_check.setCheckState(QtCore.Qt.Checked)
        vbox.addWidget(dst_grid_check)
        # dst field
        dst_field_check = QtWidgets.QCheckBox("Show field")
        dst_field_check.stateChanged.connect(self.toggle_dst_field)
        dst_field_check.setCheckState(QtCore.Qt.Checked)
        vbox.addWidget(dst_field_check)
        # set group layout
        dst_group.setLayout(vbox)
        
        show_layout.addWidget(src_group)
        show_layout.addWidget(dst_group)
        show_layout.addStretch(1)
        
        return show_widget
    
    def readCube(self, filename):
        cube = None
        cubes = iris.load(filename)
        for cb in cubes:
            if cb.var_name == 'pointData':
                cube = cb
        
        return cube
    
    def createField(self, polydata_valid, polydata_invalid, radius, scalar_bar=False):
        """Create the field."""
        actors_list = []
        
        #### valid ####
        # source
        source_valid = vtk.vtkSphereSource()
        source_valid.SetPhiResolution(self.resolution)
        source_valid.SetThetaResolution(self.resolution)
        source_valid.SetRadius(radius)
        
        # lut
        lut = vtk.vtkLookupTable()
        lut.SetNumberOfColors(1024)
        lut.SetHueRange(0.667,0.0)
        lut.SetRange(self.field_min, self.field_max)    
        lut.SetRampToLinear()
        lut.Build()
        
        # glyph mapper
        mapper_valid = vtk.vtkGlyph3DMapper()
        mapper_valid.SetInputData(polydata_valid)
        mapper_valid.SetSourceConnection(source_valid.GetOutputPort())
        mapper_valid.ClampingOff()
        mapper_valid.SetScaleFactor(1.0)
        mapper_valid.SetScaleModeToNoDataScaling()
        mapper_valid.SetLookupTable(lut)
        mapper_valid.SetScalarRange(self.field_min, self.field_max)
        
        # actor
        actor_valid = vtk.vtkActor()
        actor_valid.SetMapper(mapper_valid)
        actors_list.append(actor_valid)
        
        #### invalid ####
        if polydata_invalid is not None:
            # source
            source_invalid = vtk.vtkCubeSource()
            source_invalid.SetXLength(1.7 * radius)
            source_invalid.SetYLength(1.7 * radius)
            source_invalid.SetZLength(1.7 * radius)
            
            # glyph mapper
            mapper_invalid = vtk.vtkGlyph3DMapper()
            mapper_invalid.SetInputData(polydata_invalid)
            mapper_invalid.SetSourceConnection(source_invalid.GetOutputPort())
            mapper_invalid.ClampingOff()
            mapper_invalid.SetScaleFactor(1.0)
            mapper_invalid.SetScaleModeToNoDataScaling()
            
            # actor
            actor_invalid = vtk.vtkActor()
            actor_invalid.SetMapper(mapper_invalid)
            actor_invalid.GetProperty().SetColor(0.6, 0.6, 0.6)
            actors_list.append(actor_invalid)
        
        # create a scalar bar
        if scalar_bar:
            scalar_bar_actor = vtk.vtkScalarBarActor()
            scalar_bar_actor.SetLookupTable(lut)
            actors_list.append(scalar_bar_actor)
        
        return actors_list
    
    def createGrid(self, cube, pt, color=(1.,1.,1.), radius=1.0, show_mesh_as_surface=False):
        n0, n1 = cube.data.shape
        sg = vtk.vtkStructuredGrid()
        sg = vtk.vtkStructuredGrid()
        sg.SetDimensions(1, n0, n1)
        sg.SetPoints(pt)
        mp, ac = None, None
        # show mesh as a surface
        if show_mesh_as_surface:
            mp = vtk.vtkDataSetMapper()
            mp.SetInputData(sg)
            ac = vtk.vtkActor()
            ac.SetMapper(mp)
        # show the grid as tubes
        ed = vtk.vtkExtractEdges()
        et = vtk.vtkTubeFilter()
        em = vtk.vtkPolyDataMapper()
        ea = vtk.vtkActor()
        et.SetRadius(0.01)
        ed.SetInputData(sg)
        et.SetInputConnection(ed.GetOutputPort())
        em.SetInputConnection(et.GetOutputPort())
        ea.SetMapper(em)
        ea.GetProperty().SetColor(color)
        return [ea, ac], et, ed, em, sg, pt, mp
    
    def makePoints(self, cube, radius):
        """Create points and polydata."""
        # extract data from cube
        field = cube.data
        n0, n1 = field.shape
        numPoints = n0 * n1
        coords = cube.coords()
        lats = coords[0].points
        lons = coords[1].points
        
        # full points
        points_full = vtk.vtkPoints()
        points_full.SetNumberOfPoints(numPoints)
        
        # points for valid and invalid
        num_invalid = np.ma.count_masked(field)
        num_valid = numPoints - num_invalid
        
        # populate points_full
        k = 0
        for i1 in range(n1):
            for i0 in range(n0):
                x = radius * math.cos(lats[i0, i1] * math.pi/180.) * math.cos(lons[i0, i1] * math.pi/180.)
                y = radius * math.cos(lats[i0, i1] * math.pi/180.) * math.sin(lons[i0, i1] * math.pi/180.)
                z = radius * math.sin(lats[i0, i1] * math.pi/180.)
                points_full.SetPoint(k, x, y, z)
                k += 1
        
        # populate valid and invalid points
        if num_invalid:
            points_valid = vtk.vtkPoints()
            points_valid.SetNumberOfPoints(num_valid)
            field_valid = vtk.vtkDoubleArray()
            points_invalid = vtk.vtkPoints()
            points_invalid.SetNumberOfPoints(num_invalid)
            count_invalid = 0
            count_valid = 0
            k = 0
            min_val = None
            max_val = None
            for i1 in range(n1):
                for i0 in range(n0):
                    # position
                    x = radius * math.cos(lats[i0, i1] * math.pi/180.) * math.cos(lons[i0, i1] * math.pi/180.)
                    y = radius * math.cos(lats[i0, i1] * math.pi/180.) * math.sin(lons[i0, i1] * math.pi/180.)
                    z = radius * math.sin(lats[i0, i1] * math.pi/180.)
                    
                    if field[i0, i1] is np.ma.masked:
                        # add to invalid
                        points_invalid.SetPoint(count_invalid, x, y, z)
                        count_invalid += 1
                    
                    else:
                        # add to valid
                        points_valid.SetPoint(count_valid, x, y, z)
                        field_valid.InsertNextValue(field[i0, i1])
                        if min_val is None or field[i0, i1] < min_val:
                            min_val = field[i0, i1]
                        if max_val is None or field[i0, i1] > max_val:
                            max_val = field[i0, i1]
                        count_valid += 1
                
                k += 1
            
            if self.field_min is None or min_val < self.field_min:
                self.field_min = min_val
            if self.field_max is None or max_val > self.field_max:
                self.field_max = max_val
            
            # create polydatas
            polydata_valid = vtk.vtkPolyData()
            polydata_valid.SetPoints(points_valid)
            polydata_valid.GetPointData().SetScalars(field_valid)
            polydata_invalid = vtk.vtkPolyData()
            polydata_invalid.SetPoints(points_invalid)
        
        # case when there is no masking
        else:
            # create vtk field array
            field_full = vtk.vtkDoubleArray()
            min_val = None
            max_val = None
            for i1 in range(n1):
                for i0 in range(n0):
                    field_full.InsertNextValue(field[i0, i1])
                    if min_val is None or field[i0, i1] < min_val:
                        min_val = field[i0, i1]
                    if max_val is None or field[i0, i1] > max_val:
                        max_val = field[i0, i1]
            if self.field_min is None or min_val < self.field_min:
                self.field_min = min_val
            if self.field_max is None or max_val > self.field_max:
                self.field_max = max_val
            
            # create polydata
            polydata_valid = vtk.vtkPolyData()
            polydata_valid.SetPoints(points_full)
            polydata_valid.GetPointData().SetScalars(field_full)
            polydata_invalid = None
        
        return points_full, polydata_valid, polydata_invalid


if __name__ == "__main__":
    import sys
    import argparse
    
    parser = argparse.ArgumentParser(description='Plot grid in 3d')
    parser.add_argument('--src_file', type=str, dest='src_file', default='src.nc',
                        help='Source data file name')
    parser.add_argument('--dst_file', type=str, dest='dst_file', default='dst.nc',
                        help='Destination data file name')
    args = parser.parse_args()
    
    # start the application
    app = QtWidgets.QApplication(sys.argv)
    mainWin = MainWindow(args.src_file, args.dst_file)
    mainWin.show()
    sys.exit(app.exec_())
