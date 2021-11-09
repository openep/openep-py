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
from pyvistaqt import BackgroundPlotter
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

        self._init_data()
        self._init_ui()
        self._add_menu_bar()
        self._create_plotter_1_dock()
        self._create_plotter_2_dock()
        self._create_canvas_3_dock()
        self._create_canvas_4_dock()
        self._add_dock_widgets()

    def _init_data(self):
        """
        Set default values for variables that can later be set by the user.
        """

        # Display only the first electrogram
        self.egm_point = np.array([0], dtype=int)

        # initial bipolar voltage colourbar limits
        self._initial_lower_limit_1 = 0
        self._initial_upper_limit_1 = 2

        # initial LAT colourbar limits
        self._initial_lower_limit_2 = 0
        self._initial_upper_limit_2 = 200

    def _init_ui(self):
        """
        Initialise the GUI.

        Lots of the layout is handled automatically by QtWidgets.QMainWindow.
        """

        self.setWindowTitle("OpenEP: The open-source solution for electrophysiology data analysis")

    def _add_menu_bar(self):
        """Add a menu bar to the GUI."""

        menubar = self.menuBar()
        file_menu = menubar.addMenu("File")
        menubar.setNativeMenuBar(True)

        self.load_case_act = QtWidgets.QAction("Open File...")
        self.load_case_act.triggered.connect(self.load_data)
        file_menu.addAction(self.load_case_act)
        file_menu.addSeparator()

    def _create_plotter_1_dock(self):
        """
        Create a dockable pyvista-qt plotter for rendering 3D maps.

        This can be used for projecting scalar fields onto the 3D surface.
        Currently, only bipolar voltage is supported.
        """

        # TODO: add support for changing the scalar field

        # Plotter 1 defaults to bipolar voltage
        self.dock_1 = DockWidget("Voltage")

        plotter_1 = BackgroundPlotter(
            show=False,
            app=QtWidgets.QApplication.instance(),
            allow_quit_keypress=False,
            line_smoothing=False,
            point_smoothing=False,
            polygon_smoothing=False,
            toolbar=False,
            editor=False,
            menu_bar=False,
            title="Voltage",
            border=False,
        )
        plotter_1.background_color = 'white'
        plotter_1.setMinimumSize(QtCore.QSize(50, 50))
        self.plotter_1 = plotter_1

        # Add thresholds to plotter 1
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

        self.dock_1.setWidget(self.plotter_1)

    def _create_plotter_2_dock(self):
        """
        Create a second dockable pyvista-qt plotter for rendering 3D maps.

        This can be used for displaying alternative scalar fields to `plotter_1`.
        Currently, only local activation time is supported.
        """

        # TODO: add support for changing the scalar field

        # Plotter 2 defaults to local activation time
        self.dock_2 = DockWidget("LAT")

        plotter_2 = BackgroundPlotter(
            show=False,
            app=QtWidgets.QApplication.instance(),
            allow_quit_keypress=False,
            line_smoothing=False,
            point_smoothing=False,
            polygon_smoothing=False,
            toolbar=False,
            editor=False,
            menu_bar=False,
            title="Voltage",
            border=False,
        )
        plotter_2.background_color = 'white'
        plotter_2.setMinimumSize(QtCore.QSize(50, 50))
        self.plotter_2 = plotter_2

        # Add thresholds to plotter 2
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

        self.dock_2.setWidget(self.plotter_2)

    def _create_canvas_3_dock(self):
        """Create a dockable widget for plotting interactive electrograms with matplotlib"""

        # Canvas 3 is for the electrograms
        self.dock_3 = DockWidget("EGMs")

        self.figure_3, self.axis_3 = plt.subplots(ncols=1, nrows=1)
        self.figure_3.set_facecolor("white")
        self.canvas_3 = FigureCanvas(self.figure_3)
        self.axis_3.axis('off')  # hide them until we have data to plot
        egm_layout = QtWidgets.QFormLayout(self)

        # Add EGM selections
        self.egm_select = QtWidgets.QLineEdit("EGMs", self.canvas_3)
        self.egm_select.setStyleSheet("background-color: white; border: 1px solid lightGray;")
        self.egm_select.setGeometry(250, 0, 50, 40)
        self.egm_select.setText(str(0))
        egm_layout.addRow("EGMs", self.egm_select)

        button_egm_select = QtWidgets.QPushButton("Select EGMs (indices of points)", self.canvas_3)
        button_egm_select.setStyleSheet("background-color: lightGray")
        button_egm_select.setGeometry(0, 0, 240, 40)
        button_egm_select.clicked.connect(self.plot_electrograms)
        egm_layout.addRow(button_egm_select)

        # Add radio buttons to select bipolar, unipolar, and reference electrograms
        egm_type_layout = QtWidgets.QHBoxLayout()
        bipolar_checkbox = QtWidgets.QCheckBox("Bipolar", self.canvas_3)
        bipolar_checkbox.setStyleSheet("color: #0485d1")
        bipolar_checkbox.setGeometry(0,45, 70, 20)
        bipolar_checkbox.setChecked(True)
        bipolar_checkbox.stateChanged.connect(self.plot_electrograms)
        egm_type_layout.addWidget(bipolar_checkbox)
        
        self.bipolar_checkbox = bipolar_checkbox
        egm_layout.addRow(egm_type_layout)

        # Create toolbar, passing canvas as first parament, parent (self, the MainWindow) as second.
        toolbar = NavigationToolbar(self.canvas_3, self.dock_3)
        canvas_layout = QtWidgets.QVBoxLayout()
        canvas_layout.addWidget(self.canvas_3)
        canvas_layout.addWidget(toolbar)

        # Create a placeholder widget to hold our toolbar and canvas.
        canvas_widget = QtWidgets.QWidget()
        canvas_widget.setLayout(canvas_layout)
        canvas_widget.setStyleSheet("border-width: 0px; border: 0px; background-color:white;")

        self.dock_3.setWidget(canvas_widget)

    def _create_canvas_4_dock(self):
        """Create a dockable widget for other matplotlib plots."""

        # TODO: Add matplotlib canvas for plotting results from other analyses

        # Canvas 4 is for plotting other analyses (e.g. histogram of voltages)
        self.dock_4 = DockWidget("Analysis")

        content4 = QtWidgets.QWidget()
        content4.setStyleSheet("background-color:white;")
        content4.setMinimumSize(QtCore.QSize(50, 50))

        self.dock_4.setWidget(content4)

    def _add_dock_widgets(self):
        """Add dockable widgets to the main window."""

        self.addDockWidget(Qt.LeftDockWidgetArea, self.dock_1)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.dock_2)
        self.tabifyDockWidget(self.dock_1, self.dock_2)
        self.addDockWidget(Qt.RightDockWidgetArea, self.dock_3)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.dock_4)
        self.tabifyDockWidget(self.dock_3, self.dock_4)

        self.setDockOptions(self.GroupedDragging | self.AllowTabbedDocks | self.AllowNestedDocks)
        self.setTabPosition(Qt.AllDockWidgetAreas, QtWidgets.QTabWidget.North)

        for dock in [self.dock_1, self.dock_2, self.dock_3, self.dock_4]:
            dock.setAllowedAreas(Qt.AllDockWidgetAreas)

    def load_data(self):

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

            self.egm_times = np.arange(self.case.electric.bipolar_egm.egm.shape[1])
            self.axis_3.axis('on')  # make sure we can see the axes now
            self.axis_3.set_xlim(self.egm_times[0]-100, self.egm_times[-1]+100)
            self.axis_3.set_ylim(-0.5, 12.5)

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

    def plot_electrograms(self):

        # Get data for new set of points
        self.egm_point = np.asarray(self.egm_select.text().split(','), dtype=int)

        self.egm_bipolar_traces, self.egm_names = openep.case.get_electrograms_at_points(
            self.case,
            within_woi=False,
            buffer=0,
            indices=self.egm_point,
            egm_type="bipolar",
            return_names=True,
            return_lat=False,
        )

        # Set up axis for new plots
        ylim = self.axis_3.get_ylim()
        xlim = self.axis_3.get_xlim()
        
        self.axis_3.cla()
        self.axis_3.set_yticklabels([])
        self.axis_3.set_ylim(ylim)
        self.axis_3.set_xlim(xlim)

        # Bipolar voltage
        if self.bipolar_checkbox.isChecked():
            _, self.axis_3.axes = openep.draw.plot_electrograms(
                self.egm_times,
                self.egm_bipolar_traces,
                names=self.egm_names,
                axis=self.axis_3.axes,
            )        

        self.canvas_3.draw()



def main():

    # Create an instance of Qapplication
    app = QtWidgets.QApplication(sys.argv)

    # Create an instance of GUI
    window = OpenEpGUI()
    window.showMaximized()

    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
