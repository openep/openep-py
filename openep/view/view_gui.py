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
GUI code for OpenEp
"""
import sys

from PyQt5 import QtWidgets as qtw
from pyvistaqt import QtInteractor
import numpy as np
from copy import deepcopy
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas
)

import openep


class OpenEpGUI(qtw.QWidget):
    def __init__(self):
        super().__init__()
        self.title = "OpenEp Application"
        self.left = 10
        self.top = 10
        self.width = 840
        self.height = 480
        self.initUI()
        self.thresholds = False
        self.minval = 0
        self.maxval = 2
        self.egm_point = np.array([0], dtype=int)

    def initUI(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)

        # # LAYOUTS
        self.mainLayout = qtw.QVBoxLayout(self)
        self.menulayout = qtw.QGridLayout(self)
        self.labelLayout = qtw.QFormLayout(self)
        self.buttonLayout = qtw.QHBoxLayout(self)

        # MenuBar
        menubar = qtw.QMenuBar()
        self.menulayout.addWidget(menubar, 0, 0)
        file_menu = menubar.addMenu("File")

        self.load_case_act = qtw.QAction("Open File...")
        self.load_case_act.triggered.connect(self.load_data)
        file_menu.addAction(self.load_case_act)
        file_menu.addSeparator()

        plot_menu = qtw.QMenu("Plots", self)
        voltage_map_plot_act = qtw.QAction("3-D Voltage Map", self)
        self.voltage_map_interpolated_plot_act = qtw.QAction(
            "Interpolated Voltage Plot", self
        )
        self.historgram_plot_act = qtw.QAction("Histogram Plot", self)
        self.egm_plot_act = qtw.QAction("EGM Plot", self)

        voltage_map_plot_act.triggered.connect(self.plot_3d_map)
        plot_menu.addAction(voltage_map_plot_act)

        self.voltage_map_interpolated_plot_act.triggered.connect(
            self.plot_interpolated_voltage
        )
        plot_menu.addAction(self.voltage_map_interpolated_plot_act)

        self.historgram_plot_act.triggered.connect(self.plot_histogram)
        plot_menu.addAction(self.historgram_plot_act)

        self.egm_plot_act.triggered.connect(self.plot_egm)
        plot_menu.addAction(self.egm_plot_act)

        file_menu.addMenu(plot_menu)

        file_menu.addSeparator()

        self.quit_app_act = qtw.QAction("Quit")
        self.quit_app_act.triggered.connect(qtw.qApp.quit)
        file_menu.addAction(self.quit_app_act)

        menubar.addMenu("Edit")
        menubar.addMenu("View")
        menubar.addMenu("Help")

        # Label
        mainLabel1 = qtw.QLabel("OpenEp Tool", self)
        self.labelLayout.addWidget(mainLabel1)
        mainLabel2 = qtw.QLabel(
            "The Open Source solution for electrophysiology data analysis", self
        )
        self.labelLayout.addWidget(mainLabel2)

        self.labelLayout.addRow(self.buttonLayout)

        # Plot
        self.plotLayout = qtw.QGridLayout()
        self.plotLayout.setColumnStretch(0, 5)
        self.plotLayout.setColumnStretch(1, 5)

        self.frame = qtw.QFrame()
        self.plotter = QtInteractor(self.frame)
        self.plotLayout.addWidget(self.plotter.interactor, 0, 0)

        # # VoltThresholds limit
        self.limitLayout = qtw.QFormLayout(self)
        self.lowerlimit = qtw.QLineEdit()
        self.upperlimit = qtw.QLineEdit()
        self.lowerlimit.setText(str(0))
        self.upperlimit.setText(str(2))
        self.limitLayout.addRow("Voltage Threshold - Lower", self.lowerlimit)
        self.limitLayout.addRow("Voltage Threshold - Upper", self.upperlimit)

        button_set_thresholds = qtw.QPushButton("Set Voltage Thresholds", self)
        button_set_thresholds.setGeometry(200, 150, 100, 40)
        button_set_thresholds.clicked.connect(self.update_voltage_threshold)
        self.limitLayout.addRow(button_set_thresholds)

        # selecting the Egm points
        self.egmselect = qtw.QLineEdit()
        self.egmselect.setText(str(0))
        self.limitLayout.addRow("Egm Point", self.egmselect)

        button_egmselect = qtw.QPushButton("Set Egm Point", self)
        button_egmselect.setGeometry(200, 150, 100, 40)
        button_egmselect.clicked.connect(self.plot_egm)
        self.limitLayout.addRow(button_egmselect)

        # Nesting Layouts and setting the main layout
        self.mainLayout.addLayout(self.menulayout)
        self.mainLayout.addLayout(self.labelLayout)
        self.mainLayout.addLayout(self.plotLayout)
        self.mainLayout.addLayout(self.limitLayout)
        self.setLayout(self.mainLayout)

    def load_data(self):
        print("Please wait: Loading Data ... ")
        # Loading file from a Dialog Box
        dlg = qtw.QFileDialog()
        dlg.setFileMode(qtw.QFileDialog.AnyFile)

        if dlg.exec_():
            self.filenames = dlg.selectedFiles()
            self.ep_case = openep.load_case(self.filenames[0])
            self.mesh = self.ep_case.create_mesh()
            self.mesh1 = deepcopy(self.mesh)
            self.mesh2 = deepcopy(self.mesh)
            self.volt = self.ep_case.fields['bip']
            self.minval = 0
            self.maxval = 2
            self.free_boundaries = openep.mesh.get_free_boundaries(self.mesh)
            self.add_mesh_kws = {
                "clim": [self.minval, self.maxval],
                "scalar_bar_args": {"label_font_size": 8}
            }

            # Interpolated voltage data
            self.voltage_data = openep.case.get_voltage_electroanatomic(self.ep_case)

    def plot_3d_map(self):

        self.minval = float(self.lowerlimit.text())
        self.maxval = float(self.upperlimit.text())
        self.add_mesh_kws["clim"] = [self.minval, self.maxval]

        self.plotter = openep.draw.draw_map(
            mesh=self.mesh1,
            field=self.volt,
            plotter=self.plotter,
            free_boundaries=False,
            add_mesh_kws=self.add_mesh_kws,
        )

        self.plotter = openep.draw.draw_free_boundaries(
            self.free_boundaries,
            colour="black",
            width=5,
            plotter=self.plotter
        )

        self.plotter.reset_camera()

    def plot_interpolated_voltage(self):

        self.minval = float(self.lowerlimit.text())
        self.maxval = float(self.upperlimit.text())
        self.add_mesh_kws["clim"] = [self.minval, self.maxval]

        # # QDock Widget
        self.frame1 = qtw.QFrame()
        self.plotter1 = QtInteractor(self.frame1)
        self.dock_plot = qtw.QDockWidget("Interpolated Voltage Plot", self)
        self.dock_plot.setFloating(False)
        self.dock_plot.setWidget(self.plotter1)

        self.plotLayout.addWidget(self.dock_plot, 0, 1)

        self.plotter1 = openep.draw.draw_map(
            mesh=self.mesh2,
            field=self.voltage_data,
            plotter=self.plotter1,
            free_boundaries=False,
            add_mesh_kws=self.add_mesh_kws,
        )

        self.plotter1 = openep.draw.draw_free_boundaries(
            self.free_boundaries,
            colour="black",
            width=5,
            plotter=self.plotter1
        )

        self.plotter1.reset_camera()

    def plot_egm(self):
        self.egm_point = np.asarray(self.egmselect.text().split(','), dtype=int)

        self.fig, self.ax = plt.subplots(ncols=1, nrows=1)
        self.fig.set_facecolor("white")

        self.egm_traces, self.egm_names, self.egm_lat = openep.case.get_electrograms_at_points(
            self.ep_case, indices=self.egm_point
        )
        times = openep.case.get_woi_times(self.ep_case)
        relative_times = openep.case.get_woi_times(self.ep_case, relative=True)
        self.fig, self.ax = openep.draw.plot_electrograms(
            relative_times,
            self.egm_traces[:, times],
            names=self.egm_names,
        )
        self.canvas = FigureCanvas(self.fig)

        self.dock_plot1 = qtw.QDockWidget("EGM Plot", self)
        self.dock_plot1.setFloating(False)
        self.dock_plot1.setWidget(self.canvas)
        self.plotLayout.addWidget(self.dock_plot1, 1, 1)

    def plot_histogram(self):
        self.plotter3 = QtInteractor(self.frame)
        self.dock_plot2 = qtw.QDockWidget("Histogram Plot", self)
        self.dock_plot2.setFloating(False)
        self.dock_plot2.setWidget(self.plotter3)

        self.plotLayout.addWidget(self.dock_plot2, 1, 0)
        pass

    def update_voltage_threshold(self):
        self.plot_3d_map()
        self.plot_interpolated_voltage()


def main():
    # Create an instance of Qapplication
    app = qtw.QApplication(sys.argv)
    # Create an instance of GUI
    window = OpenEpGUI()
    window.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
