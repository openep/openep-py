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
import pathlib

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import Qt
import numpy as np

import openep
import openep.view.custom_widgets
import openep.view.plotters
import openep.view.canvases
import openep.view.images


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
        self._disable_dock_widgets()

    def _init_data(self):
        """
        Set default values for variables that can later be set by the user.
        """

        # Initially, display first electrogram only
        self.egm_points = np.array([0], dtype=int)

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

        self.menubar = self.menuBar()
        file_menu = self.menubar.addMenu("File")
        self.menubar.setNativeMenuBar(True)

        self.load_case_act = QtWidgets.QAction("Open File...")
        self.load_case_act.triggered.connect(self._load_data)
        file_menu.addAction(self.load_case_act)
        file_menu.addSeparator()

    def _create_plotter_1_dock(self):
        """
        Create a dockable pyvista-qt plotter for rendering 3D maps.

        This can be used for projecting the voltage the 3D surface.
        The user can select to use the clinical values or values interpolated
        from the bipolar or unipolar electrograms.
        """

        # Plotter 1 is used for the bipolar voltage
        self.dock_1 = openep.view.custom_widgets.CustomDockWidget("Voltage")
        self.plotter_1 = openep.view.plotters.create_plotter()

        # Add thresholds to plotter 1 colourbar
        self.plotter_1.plotterLayout = QtWidgets.QFormLayout(self)
        self.lower_limit_1, self.upper_limit_1, button_set_thresholds_1 = \
            openep.view.plotters.create_colourbar_widgets(
                self.plotter_1,
                lower=self._initial_lower_limit_1,
                upper=self._initial_upper_limit_1,
                )
        button_set_thresholds_1.clicked.connect(self.update_colourbar_limits_1)

        self.plotter_1.plotterLayout.addRow(self.lower_limit_1)
        self.plotter_1.plotterLayout.addRow(self.upper_limit_1)
        self.plotter_1.plotterLayout.addRow(button_set_thresholds_1)

        # Add radio buttons to select clinical or interpolated bipolar voltages
        map_type_layout = QtWidgets.QHBoxLayout()
        radio_buttons = openep.view.plotters.create_map_type_widgets(self.plotter_1)
        self.plotter_1_clinical_radio = radio_buttons[0]
        self.plotter_1_openep_bipolar_radio = radio_buttons[1]
        self.plotter_1_openep_unipolar_radio = radio_buttons[2]
        self.plotter_1_radio_checked = self.plotter_1_clinical_radio  # default checked button

        self.plotter_1_clinical_radio.toggled.connect(
            lambda: self.set_plotter_1_button_state(self.plotter_1_clinical_radio)
        )

        self.plotter_1_openep_bipolar_radio.toggled.connect(
            lambda: self.set_plotter_1_button_state(self.plotter_1_openep_bipolar_radio)
        )
        self.plotter_1_openep_unipolar_radio.toggled.connect(
            lambda: self.set_plotter_1_button_state(self.plotter_1_openep_unipolar_radio)
        )

        map_type_layout.addWidget(self.plotter_1_clinical_radio)
        map_type_layout.addWidget(self.plotter_1_openep_bipolar_radio)
        map_type_layout.addWidget(self.plotter_1_openep_unipolar_radio)
        self.plotter_1.plotterLayout.addRow(map_type_layout)

        self.dock_1.setWidget(self.plotter_1)

    def _create_plotter_2_dock(self):
        """
        Create a dockable pyvista-qt plotter for rendering 3D maps.

        This can be used for projecting the local activation time the 3D surface.
        Currently, only the values obtained from the clinical mapping system can
        be used (interpolation based on electrograms is not yet supported.)
        """

        # Plotter 2 is used for the local activation time
        self.dock_2 = openep.view.custom_widgets.CustomDockWidget("LAT")
        self.plotter_2 = openep.view.plotters.create_plotter()

        # Add thresholds to plotter 1
        self.plotter_2.plotterLayout = QtWidgets.QFormLayout(self)
        self.lower_limit_2, self.upper_limit_2, button_set_thresholds_2 = \
            openep.view.plotters.create_colourbar_widgets(
                self.plotter_2,
                lower=self._initial_lower_limit_2,
                upper=self._initial_upper_limit_2,
                )
        button_set_thresholds_2.clicked.connect(self.update_colourbar_limits_2)

        self.plotter_2.plotterLayout.addRow(self.lower_limit_2)
        self.plotter_2.plotterLayout.addRow(self.upper_limit_2)
        self.plotter_2.plotterLayout.addRow(button_set_thresholds_2)

        self.dock_2.setWidget(self.plotter_2)

    def _create_canvas_3_dock(self):
        """
        Create a dockable widget for plotting interactive electrograms with matplotlib.

        The user can select which points the electrograms will be plotted for, as well
        as the type(s) of electrograms to plot: reference, bipolar, unipolar A,
        unipolar B.

        The user can also change the window of interest with a RangeSlider,
        and setting the window of interest with a push button will re-interpolate
        the electrogram data onto the surface.
        """

        self.dock_3 = openep.view.custom_widgets.CustomDockWidget("EGMs")
        self.canvas_3, self.figure_3, self.axis_3 = openep.view.canvases.create_canvas()
        egm_layout = QtWidgets.QFormLayout(self)

        # Add EGM selections
        self.egm_select, button_egm_select = openep.view.canvases.add_egm_select(
            canvas=self.canvas_3,
        )
        button_egm_select.clicked.connect(self.update_electrograms)
        egm_layout.addRow("EGMs", self.egm_select)
        egm_layout.addRow(button_egm_select)

        # Add radio buttons to select bipolar, unipolar, and reference electrograms
        egm_type_layout = QtWidgets.QHBoxLayout()
        self.reference_checkbox, self.bipolar_checkbox, self.unipolar_A_checkbox, self.unipolar_B_checkbox = \
            openep.view.canvases.add_egm_type_widgets(
                canvas=self.canvas_3,
            )

        for box in [self.reference_checkbox, self.bipolar_checkbox,
                    self.unipolar_A_checkbox, self.unipolar_B_checkbox]:

            box.stateChanged.connect(self.plot_electrograms)
            egm_type_layout.addWidget(box)

        egm_layout.addRow(egm_type_layout)

        # Create toolbar
        toolbar = openep.view.canvases.add_toolbar(
            canvas=self.canvas_3,
            parent=self.dock_3,
        )

        # Create a placeholder widget to hold our toolbar and canvas.
        canvas_widget = openep.view.canvases.create_canvas_widget(
            canvas=self.canvas_3,
            toolbar=toolbar,
        )

        self.dock_3.setWidget(canvas_widget)

    def _create_canvas_4_dock(self):
        """
        Create a dockable widget for other matplotlib plots.

        For example, plotting a histogram of the surface area occupied by
        a range of bipolar voltages.
        """

        self.dock_4 = openep.view.custom_widgets.CustomDockWidget("Analysis")
        self.canvas_4, self.figure_4, self.axis_4 = openep.view.canvases.create_canvas()
        toolbar = openep.view.canvases.add_toolbar(
            canvas=self.canvas_4,
            parent=self.dock_4,
            keep_actions=['Save'],
        )

        # Create a placeholder widget to hold our toolbar and canvas.
        canvas_widget = openep.view.canvases.create_canvas_widget(
            canvas=self.canvas_4,
            toolbar=toolbar,
        )

        self.dock_4.setWidget(canvas_widget)

    def _add_dock_widgets(self):
        """
        Add dockable widgets to the main window.

        The two BackgroundPlotters are tabified, as are the two MPL canvases.
        """

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

    def _disable_dock_widgets(self):
        """
        Ignore all key presses in the dock widgets.

        This is required if there is no file loaded, otherwise
        the GUI will crash.
        """

        self.dock_1.setEnabled(False)
        self.dock_2.setEnabled(False)
        self.dock_3.setEnabled(False)
        self.dock_4.setEnabled(False)

    def _enable_dock_widgets(self, exclude=None):
        """
        Enable all dock widgets.

        Once a file has been loaded it is save to enable to dock widgets.

        Args:
            exclude (list): a list of dock widgets to exclude from enabling.
                May be useful if a different file format is loaded and it
                does not contain data required by some dock widgets.
        """
        exclude = exclude if exclude is not None else []

        for dock in [self.dock_1, self.dock_2, self.dock_3, self.dock_4]:
            if dock not in exclude:
                dock.setEnabled(True)

    def _load_data(self):
        """Load an OpenEP case."""

        dialogue = QtWidgets.QFileDialog()
        dialogue.setWindowTitle('Load an OpenEP dataset')
        dialogue.setDirectory(QtCore.QDir.currentPath())
        dialogue.setFileMode(QtWidgets.QFileDialog.ExistingFiles)
        dialogue.selectNameFilter("MATLAB files (*.mat);;openCARP files (*.pts *.elem *.dat)")
        dialogue.setNameFilters(["MATLAB files (*.mat)", "openCARP files (*.pts *.elem *.dat)"])

        if dialogue.exec_():

            filenames = dialogue.selectedFiles()
            if (len(filenames) == 0) and (pathlib.Path(filenames[0]).suffix == ".mat"):
                self._load_case_and_initialise(filenames[0])
            elif len(filenames) == 3:
                self._load_openCARP_and_initialise(filenames)
            else:

                error = QtWidgets.QMessageBox()
                error.setIcon(QtWidgets.QMessageBox.Critical)
                error.setText("File selection error")
                error.setInformativeText(
                    "Please select either a single MATLAB file (*.mat), or a single "
                    "set of openCARP files (one each of *.pts, *.elem, *.dat)."
                )
                error.setWindowTitle("Error")
                error.exec_()

    def _load_openCARP_and_initialise(self, filenames):
        """Load data from an openCARP simulation.

        """

        # TODO: create a new dock widget for the simulation data
        # Allow user to draw either bipolar or unipolar voltages

        extensions = [pathlib.Path(f).suffix for f in filenames]
        if (".pts" in extensions) and (".elem" in extensions) and (".dat" in extensions):
            self.carp = openep.load_openCARP(
                points=filenames[extensions.index(".pts")],
                indices=filenames[extensions.index(".elem")],
                unipolar_egm=filenames[extensions.index(".dat")],
            )

    def _load_case_and_initialise(self, filename):
        """Load an OpenEP case object.

        Create separate meshes for the two BackgroundPlotters.
        Interpolate data ready to be plotted.
        """

        self.case = openep.load_case(filename)
        self.mesh_1 = self.case.create_mesh()
        self.mesh_2 = self.case.create_mesh()
        self.free_boundaries = openep.mesh.get_free_boundaries(self.mesh_1)

        self.add_mesh_1_kws = {
            "clim": [self._initial_lower_limit_1, self._initial_upper_limit_1],
            "scalar_bar_args": {
                "title": "Voltage (mV)",
                "title_font_size": 8,
                "label_font_size": 8,
                "color": "#363737",
                "vertical": False,
                "width": 0.4,
                "height": 0.04,
                "position_x": 0.025,
                "interactive": False,
                "below_label": " ",
                "above_label": " ",
            }
        }
        self.add_mesh_2_kws = {
            "clim": [self._initial_lower_limit_2, self._initial_upper_limit_2],
            "scalar_bar_args": {
                "title": "Local activation time (ms)",
                "title_font_size": 8,
                "label_font_size": 8,
                "color": "#363737",
                "vertical": False,
                "width": 0.4,
                "height": 0.04,
                "position_x": 0.025,
                "interactive": False,
                "below_label": " ",
                "above_label": " ",
            }
        }

        # Interpolate some fields (e.g. bipolar volage)
        self.interpolate_fields()

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
        self.axis_3.set_ylim(-1, 13)

        # We can't create the slider earlier because we need to know the min and max egm times
        self.slider = openep.view.canvases.add_woi_range_slider(
            figure=self.figure_3,
            axis=[0.1535, 0.05, 0.7185, 0.01],  # location of the RangeSlider on the figure
            valmin=self.egm_times[0],
            valmax=self.egm_times[-1],
        )
        self.initialise_woi_slider_limits()
        self.slider.on_changed(self.update_woi_slider_limits)

        self.set_woi_button = openep.view.canvases.add_woi_set_button(
            figure=self.figure_3,
            axis=[0.782, 0.02, 0.09, 0.025],
        )
        self.set_woi_button.on_clicked(self.interpolate_fields_and_draw)

        # Draw default electrograms - bipolar and unipolar egms for point at index 0
        self.update_electrograms()

        # Un-freeze the GUI
        self._enable_dock_widgets()

    def initialise_woi_slider_limits(self):
        """Set the limits to be the window of interest."""

        start_woi, stop_woi = self.case.electric.annotations.window_of_interest[0]
        start_woi += self.case.electric.annotations.reference_activation_time[0]
        stop_woi += self.case.electric.annotations.reference_activation_time[0]

        self.slider.set_val([start_woi, stop_woi])

        self.slider_lower_limit = self.axis_3.axvline(
            start_woi,
            color="grey",
            linestyle='--',
            linewidth=0.8,
            alpha=0.6,
        )

        self.slider_upper_limit = self.axis_3.axvline(
            stop_woi,
            color="grey",
            linestyle='--',
            linewidth=0.8,
            alpha=0.6,
        )

    def update_woi_slider_limits(self, val):
        """
        Take the min and max values from the RangeSlider widget.
        Use this to set the window of interest and to change the x location
        of the two axvlines drawn on the EGM canvas.
        """

        # from https://matplotlib.org/devdocs/gallery/widgets/range_slider.html
        start_woi, stop_woi = val
        self.slider_lower_limit.set_xdata([start_woi, start_woi])
        self.slider_upper_limit.set_xdata([stop_woi, stop_woi])

        reference_annotation = self.case.electric.annotations.reference_activation_time[0]
        self.case.electric.annotations.window_of_interest[:, 0] = start_woi - reference_annotation
        self.case.electric.annotations.window_of_interest[:, 1] = stop_woi - reference_annotation

    def update_colourbar_limits_1(self):

        lower_limit = float(self.lower_limit_1.text())
        upper_limit = float(self.upper_limit_1.text())
        self.add_mesh_1_kws["clim"] = [lower_limit, upper_limit]

        # Now we need to redraw the map using the correct scalars
        # TODO: this should be replaced to use pyvista's set active scalars
        self.set_plotter_1_button_state(self.plotter_1_radio_checked)

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

    def interpolate_fields_and_draw(self, event=None):
        """
        Interpolate EGM data onto the surface and draw a map if necessary.

        The event argument is ignored. It is there because when
        self.set_woi_button.on_clicked (mpl.widgets.Button) is pressed,
        matplotlib passes an event to the called function (i.e. this one).
        """

        self.interpolate_fields()

        if self.plotter_1_openep_bipolar_radio.isChecked():
            self.draw_map(
                mesh=self.mesh_1,
                plotter=self.plotter_1,
                data=self.interpolated_fields['bipolar_voltage'],
                add_mesh_kws=self.add_mesh_1_kws
            )
        if self.plotter_1_openep_unipolar_radio.isChecked():
            self.draw_map(
                mesh=self.mesh_1,
                plotter=self.plotter_1,
                data=self.interpolated_fields['unipolar_voltage'],
                add_mesh_kws=self.add_mesh_1_kws
            )

    def interpolate_fields(self):
        """
        Interpolate EGM data onto the surface.

        We use a buffer of zero as the window of interest if determined by the user,
        via self.slider (mpl.widgets.RangeSlider) and self.set_woi_button (mpl.widgets.Button).
        """

        bipolar_voltage = openep.case.interpolate_voltage_onto_surface(
            self.case,
            max_distance=None,
            buffer=0,
        )
        unipolar_voltage = openep.case.interpolate_voltage_onto_surface(
            self.case,
            max_distance=None,
            buffer=0,
            bipolar=False,
        )

        self.interpolated_fields = {
            "bipolar_voltage": bipolar_voltage,
            "unipolar_voltage": unipolar_voltage,
        }

    def set_plotter_1_button_state(self, button):
        """
        Draw a bipolar voltage map.
        The type of map drawn (clinical or interpolated from electrograms) depends
        on which button was pressed.
        """

        if (button.text() == "Clinical") and button.isChecked():
            self.draw_map(
                mesh=self.mesh_1,
                plotter=self.plotter_1,
                data=self.case.fields.bipolar_voltage,
                add_mesh_kws=self.add_mesh_1_kws
            )
            self.plotter_1_radio_checked = self.plotter_1_clinical_radio

        elif (button.text() == "OpenEP - Bipolar") and button.isChecked():
            self.draw_map(
                mesh=self.mesh_1,
                plotter=self.plotter_1,
                data=self.interpolated_fields['bipolar_voltage'],
                add_mesh_kws=self.add_mesh_1_kws
            )
            self.plotter_1_radio_checked = self.plotter_1_openep_bipolar_radio

        elif (button.text() == "OpenEP - Unipolar") and button.isChecked():
            self.draw_map(
                mesh=self.mesh_1,
                plotter=self.plotter_1,
                data=self.interpolated_fields['unipolar_voltage'],
                add_mesh_kws=self.add_mesh_1_kws
            )
            self.plotter_1_radio_checked = self.plotter_1_openep_unipolar_radio

    def draw_map(self, mesh, plotter, data, add_mesh_kws):
        """
        Project scalar values onto the surface of the mesh.

        Also draw the free boundaries (anatomical structures) as black lines.
        """

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

    def update_electrograms(self):
        """
        Extract electrograms at specified indices and re-plot.

        The points are specified by the user and set via use of the
        self.egm_select (QLineEdit) and button_egm_select (QPushButton) widgets.
        """

        # Get data for new set of points
        self.egm_points = np.asarray(self.egm_select.text().split(','), dtype=int)

        self.extract_electrograms()
        self.plot_electrograms()

    def extract_electrograms(self):
        """Extract electrograms at specified indices"""

        self.egm_reference_traces, self.egm_names = openep.case.get_electrograms_at_points(
            self.case,
            within_woi=False,
            buffer=0,
            indices=self.egm_points,
            egm_type="reference",
            return_names=True,
            return_lat=False,
        )

        self.egm_bipolar_traces = openep.case.get_electrograms_at_points(
            self.case,
            within_woi=False,
            buffer=0,
            indices=self.egm_points,
            egm_type="bipolar",
            return_names=False,
            return_lat=False,
        )

        unipolar_traces = openep.case.get_electrograms_at_points(
            self.case,
            within_woi=False,
            buffer=0,
            indices=self.egm_points,
            egm_type="unipolar",
            return_names=False,
            return_lat=False,
        )
        self.egm_unipolar_A_traces = unipolar_traces[:, :, 0]
        self.egm_unipolar_B_traces = unipolar_traces[:, :, 1]

    def plot_electrograms(self):
        """
        Plot electrograms for the currently-selected set of points.

        Electrograms must first have been extracted using update_electrograms.
        Here we will plot the reference, bipolar, unipolar A, unipolar B
        electrograms for each point.
        """

        # Set up axis for new plots
        ylim = self.axis_3.get_ylim()
        xlim = self.axis_3.get_xlim()

        self.axis_3.cla()
        self.axis_3.set_yticklabels([])
        self.axis_3.set_ylim(ylim)
        self.axis_3.set_xlim(xlim)

        # Reference voltage
        if self.reference_checkbox.isChecked():

            _, self.axis_3.axes = openep.draw.plot_electrograms(
                self.egm_times,
                self.egm_reference_traces,
                axis=self.axis_3.axes,
                colour="xkcd:scarlet",
                y_separation=2,
            )

        # Bipolar voltage
        if self.bipolar_checkbox.isChecked():

            _, self.axis_3.axes = openep.draw.plot_electrograms(
                self.egm_times,
                self.egm_bipolar_traces,
                axis=self.axis_3.axes,
                colour="xkcd:cerulean",
                y_separation=2,
            )

        # Unipolar A voltage
        if self.unipolar_A_checkbox.isChecked():

            _, self.axis_3.axes = openep.draw.plot_electrograms(
                self.egm_times,
                self.egm_unipolar_A_traces,
                axis=self.axis_3.axes,
                colour="xkcd:tree green",
                y_start=1,
                y_separation=2,
            )

        # Unipolar B voltage
        if self.unipolar_B_checkbox.isChecked():

            _, self.axis_3.axes = openep.draw.plot_electrograms(
                self.egm_times,
                self.egm_unipolar_B_traces,
                axis=self.axis_3.axes,
                colour="xkcd:pumpkin",
                y_start=1,
                y_separation=2,
            )

        # Add labels if necessary
        yticks = []
        yticklabels = []
        separations = np.arange(self.egm_points.size) * 2

        if self.reference_checkbox.isChecked() or self.bipolar_checkbox.isChecked():
            yticks.extend(separations)
            yticklabels.extend(self.egm_names)

        # Unipolar A and B are shifted above the bipolar and reference electrograms for clarity
        if self.unipolar_A_checkbox.isChecked() or self.unipolar_B_checkbox.isChecked():
            yticks.extend(separations + 1)
            yticklabels.extend(self.egm_names)

        self.axis_3.set_yticks(yticks)
        self.axis_3.set_yticklabels(yticklabels)

        # draw vertical lines at the window of interest
        self.initialise_woi_slider_limits()

        self.canvas_3.draw()


def main():

    # Create an instance of Qapplication
    app = QtWidgets.QApplication(sys.argv)
    app.setWindowIcon(QtGui.QIcon(openep.view.images.LOGO))

    # Create an instance of GUI
    window = OpenEpGUI()
    window.showMaximized()

    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
