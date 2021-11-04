# OpenEP
# Copyright (c) 2021 OpenEP Collaborators
#
# This file is part of OpenEP.
#
# OpenEP is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# OpenEP is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program (LICENSE.txt).  If not, see <http://www.gnu.org/licenses/>


"""
A GUI for OpenEP-Py.
"""
import sys

from PyQt5 import QtCore, QtWidgets
from PyQt5.QtCore import Qt
from pyvistaqt import QtInteractor
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar
)

import numpy as np
import matplotlib.pyplot as plt

import openep


class DockWidget(QtWidgets.QDockWidget):
    def __init__(self, title: str):
        super().__init__(title)
        self.setTitleBarWidget(QtWidgets.QWidget())
        self.dockLocationChanged.connect(self.on_dockLocationChanged)

    def on_dockLocationChanged(self):
        main: QtWidgets.QMainWindow = self.parent()
        all_dock_widgets = main.findChildren(QtWidgets.QDockWidget)

        for dock_widget in all_dock_widgets:
            sibling_tabs = main.tabifiedDockWidgets(dock_widget)
            # If you pull a tab out of a group the other tabs still see it as a sibling while dragging...
            sibling_tabs = [s for s in sibling_tabs if not s.isFloating()]

            if len(sibling_tabs) != 0:
                # Hide title bar
                dock_widget.setTitleBarWidget(QtWidgets.QWidget())
            else:
                # Re-enable title bar
                dock_widget.setTitleBarWidget(None)

    def minimumSizeHint(self) -> QtCore.QSize:
        return QtCore.QSize(100, 100)


class MplCanvas(FigureCanvas):

    def __init__(self, nrows=1, ncols=1):
        figure, self.axes = plt.subplots(nrows=nrows, ncols=ncols)
        super(MplCanvas, self).__init__(figure)


class OpenEpGUI(QtWidgets.QMainWindow):

    def __init__(self):

        super().__init__()
        self._initUI()
        self._add_menu_bar()
        
        self.egm_point = np.array([0], dtype=int)

        # initial bipolar voltage colourbar limits
        self._initial_lower_limit_1 = 0
        self._initial_upper_limit_1 = 2
        
        # initial LAT colourbar limits
        self._initial_lower_limit_2 = 0
        self._initial_upper_limit_2 = 200

        self._add_dock_widgets()

    def _initUI(self):

        self.setWindowTitle("OpenEP: The open-source solution for electrophysiology data analysis")
        self.mainLayout = QtWidgets.QVBoxLayout()
        self.menulayout = QtWidgets.QGridLayout()
        self.buttonLayout = QtWidgets.QHBoxLayout()
    
    def _add_menu_bar(self):

        menubar = QtWidgets.QMenuBar()
        self.menulayout.addWidget(menubar, 0, 0)
        file_menu = menubar.addMenu("File")
        
        self.load_case_act = QtWidgets.QAction("Open File...")
        self.load_case_act.triggered.connect(self.load_data)
        file_menu.addAction(self.load_case_act)
        file_menu.addSeparator()
        
    def _add_dock_widgets(self):

        dock_1 = DockWidget("Voltage")
        dock_2 = DockWidget("LAT")
        dock_3 = DockWidget("EGMs")
        dock_4 = DockWidget("Analysis")

        # Plotter 1 (defaults to bipolar voltage)
        frame_1 = QtWidgets.QFrame()
        plotter_1 = QtInteractor(frame_1)
        plotter_1.background_color = 'white'
        plotter_1.setMinimumSize(QtCore.QSize(50, 50))
        self.plotter_1 = plotter_1
        
        # Add thresholds to plotter 1
        self.plotter_1.limitLayout = QtWidgets.QFormLayout(self)
        
        self.lower_limit_1 = QtWidgets.QLineEdit("Lower threshold", self.plotter_1)
        self.lower_limit_1.setStyleSheet("background-color: white; border: 1px solid lightGray;")
        self.lower_limit_1.setGeometry(200, 10, 50, 40)
        self.lower_limit_1.setText(str(self._initial_lower_limit_1))
        self.plotter_1.limitLayout.addRow("Lower threshold", self.lower_limit_1)
        
        self.upper_limit_1 = QtWidgets.QLineEdit("Upper threshold", self.plotter_1)
        self.upper_limit_1.setStyleSheet("background-color: white; border: 1px solid lightGray;")
        self.upper_limit_1.setGeometry(260, 10, 50, 40)
        self.upper_limit_1.setText(str(self._initial_upper_limit_1))
        self.plotter_1.limitLayout.addRow("Upper threshold", self.upper_limit_1)

        button_set_thresholds_1 = QtWidgets.QPushButton("Set colourbar limits:", self.plotter_1)
        button_set_thresholds_1.setStyleSheet("background-color: lightGray")
        button_set_thresholds_1.setGeometry(10, 10, 180, 40)
        button_set_thresholds_1.clicked.connect(self.update_colourbar_limits_1)
        self.plotter_1.limitLayout.addRow(button_set_thresholds_1)
        
        # Plotter 2 (defaults to LAT)
        frame_2 = QtWidgets.QFrame()
        plotter_2= QtInteractor(frame_2)
        plotter_2.background_color = 'white'
        plotter_2.setMinimumSize(QtCore.QSize(50, 50))
        self.plotter_2 = plotter_2

        # Add thresholds to plotter 2
        self.plotter_2.limitLayout = QtWidgets.QFormLayout(self)        
        
        self.lower_limit_2 = QtWidgets.QLineEdit("Lower threshold", self.plotter_2)
        self.lower_limit_2.setStyleSheet("background-color: white; border: 1px solid lightGray;")
        self.lower_limit_2.setGeometry(200, 10, 50, 40)
        self.lower_limit_2.setText(str(self._initial_lower_limit_2))
        self.plotter_2.limitLayout.addRow("Lower threshold", self.lower_limit_2)
        
        self.upper_limit_2 = QtWidgets.QLineEdit("Upper threshold", self.plotter_2)
        self.upper_limit_2.setStyleSheet("background-color: white; border: 1px solid lightGray;")
        self.upper_limit_2.setGeometry(260, 10, 50, 40)
        self.upper_limit_2.setText(str(self._initial_upper_limit_2))
        self.plotter_2.limitLayout.addRow("Upper threshold", self.upper_limit_2)

        button_set_thresholds_2 = QtWidgets.QPushButton("Set colourbar limits:", self.plotter_2)
        button_set_thresholds_2.setStyleSheet("background-color: lightGray")
        button_set_thresholds_2.setGeometry(10, 10, 180, 40)
        button_set_thresholds_2.clicked.connect(self.update_colourbar_limits_2)
        self.plotter_2.limitLayout.addRow(button_set_thresholds_2)

        # MPL canvases
        # Canvas 1 (EGMs)
        self.figure_3, self.axis_3 = plt.subplots(ncols=1, nrows=1)
        self.figure_3.set_facecolor("white")
        self.canvas_3 = FigureCanvas(self.figure_3)
        self.axis_3.axis('off')  # hide them until we have data to plot

        # Add EGM selections
        egm_layout = QtWidgets.QFormLayout(self)
        
        self.egm_select = QtWidgets.QLineEdit("EGMs", self.canvas_3)
        self.egm_select.setStyleSheet("background-color: white; border: 1px solid lightGray;")
        self.egm_select.setGeometry(260, 10, 50, 40)
        self.egm_select.setText(str(0))
        egm_layout.addRow("EGMs", self.egm_select)

        button_egm_select = QtWidgets.QPushButton("Select EGMs (indices of points)", self.canvas_3)
        button_egm_select.setStyleSheet("background-color: lightGray")
        button_egm_select.setGeometry(10, 10, 240, 40)
        button_egm_select.clicked.connect(self.plot_electrograms)
        egm_layout.addRow(button_egm_select)
        
        # Create toolbar, passing canvas as first parament, parent (self, the MainWindow) as second.
        toolbar = NavigationToolbar(self.canvas_3, dock_3)
        canvas_layout = QtWidgets.QVBoxLayout()
        canvas_layout.addWidget(self.canvas_3)
        canvas_layout.addWidget(toolbar)

        # Create a placeholder widget to hold our toolbar and canvas.
        canvas_widget = QtWidgets.QWidget()
        canvas_widget.setLayout(canvas_layout)
        canvas_widget.setStyleSheet("border-width: 0px; border: 0px; background-color:white;")

        # Other analysis (e.g. histogram of voltages)
        content4 = QtWidgets.QWidget()
        content4.setStyleSheet("background-color:white;")
        content4.setMinimumSize(QtCore.QSize(50, 50))

        dock_1.setWidget(self.plotter_1)
        dock_2.setWidget(self.plotter_2)
        dock_3.setWidget(canvas_widget)
        dock_4.setWidget(content4)
        
        self.addDockWidget(Qt.LeftDockWidgetArea, dock_1)
        self.addDockWidget(Qt.LeftDockWidgetArea, dock_2)
        self.tabifyDockWidget(dock_1, dock_2)
        self.addDockWidget(Qt.RightDockWidgetArea, dock_3)
        self.addDockWidget(Qt.LeftDockWidgetArea, dock_4)
        self.tabifyDockWidget(dock_3, dock_4)

        self.setDockOptions(self.GroupedDragging | self.AllowTabbedDocks | self.AllowNestedDocks)
        self.setTabPosition(Qt.AllDockWidgetAreas, QtWidgets.QTabWidget.North)

        for dock in [dock_1, dock_2, dock_3, dock_4]:
            dock.setAllowedAreas(Qt.AllDockWidgetAreas)
        
    def load_data(self):
        print("Please wait: Loading Data ... ")

        dialogue = QtWidgets.QFileDialog()
        dialogue.setWindowTitle('Load an OpenEP dataset')
        dialogue.setDirectory(QtCore.QDir.currentPath())
        dialogue.setFileMode(QtWidgets.QFileDialog.ExistingFile)
        dialogue.setNameFilter("MATLAB files (*.mat)")

        if dialogue.exec_():

            filename = dialogue.selectedFiles()[0]
            self.case = openep.load_case(filename)
            self.mesh_1 = self.case.create_mesh()
            self.mesh_2 = self.case.create_mesh()
            self.free_boundaries = openep.mesh.get_free_boundaries(self.mesh_1)

            self.add_mesh_1_kws = {
                "clim": [self._initial_lower_limit_1, self._initial_upper_limit_1],
                "scalar_bar_args": {"label_font_size": 8}
            }
            self.add_mesh_2_kws = {
                "clim": [self._initial_lower_limit_2, self._initial_upper_limit_2],
                "scalar_bar_args": {"label_font_size": 8}
            }

            # Interpolate some fields (e.g. bipolar volage)
            interpolated_voltage = openep.case.interpolate_voltage_onto_surface(
                self.case,
                max_distance=None,
            )
            self.interpolated_fields = {"bipolar_voltage": interpolated_voltage}
            
            # Draw default maps
            self.draw_map(
                mesh=self.mesh_1,
                plotter=self.plotter_1,
                data=self.case.fields.bipolar_voltage,
                add_mesh_kws=self.add_mesh_1_kws,
            )
            self.draw_map(
                mesh=self.mesh_2,
                plotter=self.plotter_2,
                data=self.case.fields.local_activation_time,
                add_mesh_kws=self.add_mesh_2_kws,
            )

            self.axis_3.axis('on')  # make sure we can see the axes now
            self.plot_electrograms()
    
    def update_colourbar_limits_1(self):
        
        lower_limit = float(self.lower_limit_1.text())
        upper_limit = float(self.upper_limit_1.text())
        self.add_mesh_1_kws["clim"] = [lower_limit, upper_limit]
        self.draw_map(
            mesh=self.mesh_1,
            plotter=self.plotter_1,
            data=self.case.fields.bipolar_voltage,
            add_mesh_kws=self.add_mesh_1_kws,
        )

    def update_colourbar_limits_2(self):
        
        lower_limit = float(self.lower_limit_2.text())
        upper_limit = float(self.upper_limit_2.text())
        self.add_mesh_2_kws["clim"] = [lower_limit, upper_limit]
        self.draw_map(
            mesh=self.mesh_2,
            plotter=self.plotter_2,
            data=self.case.fields.local_activation_time,
            add_mesh_kws=self.add_mesh_2_kws
        )

    def draw_map(self, mesh, plotter, data, add_mesh_kws):

        plotter = openep.draw.draw_map(
            mesh=mesh,
            field=data,
            plotter=plotter,
            free_boundaries=False,
            add_mesh_kws=add_mesh_kws,
        )

        plotter = openep.draw.draw_free_boundaries(
            self.free_boundaries,
            colour="black",
            width=5,
            plotter=plotter,
        )

        plotter.reset_camera()

    def plot_electrograms(self):

        self.egm_point = np.asarray(self.egm_select.text().split(','), dtype=int)

        self.egm_traces, self.egm_names, self.egm_lat = openep.case.get_electrograms_at_points(
            self.case, indices=self.egm_point
        )
        times = openep.case.get_woi_times(self.case)
        relative_times = openep.case.get_woi_times(self.case, relative=True)
        
        _, self.axis_3.axes = openep.draw.plot_electrograms(
            relative_times,
            self.egm_traces[:, times],
            names=self.egm_names,
            axis=self.axis_3.axes,
        )
        self.canvas_3.draw()


if __name__ == "__main__":

    app = QtWidgets.QApplication(sys.argv)
    main = OpenEpGUI()
    main.showMaximized()
    sys.exit(app.exec_())
