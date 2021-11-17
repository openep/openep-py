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
from typing import Union
from attr import attrs

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import Qt
import pyvista
import numpy as np

import openep
import openep.view.custom_widgets
import openep.view.plotters
import openep.view.canvases
import openep.view.images


class OpenEpGUI(QtWidgets.QMainWindow):

    def __init__(self):

        super().__init__()

        self._init_systems()
        self._init_ui()
        self._add_menu_bar()

        self.setDockOptions(self.GroupedDragging | self.AllowTabbedDocks | self.AllowNestedDocks)
        self.setTabPosition(Qt.AllDockWidgetAreas, QtWidgets.QTabWidget.North)

    def _init_systems(self):
        """Containers and variables for keeping track of the systems loaded into the GUI."""

        # We need to keep track of all cases currently loaded (both OpenEP and openCARP)
        # Each system will keep track of its own docks/plotters/meshes
        self.systems = []
        self._system_counter = 0

        # We will use this to determine which egms to plot
        # All windows belonging to the main case will be given red borders
        self._active_system = None

    def _init_ui(self):
        """
        Initialise the GUI.

        Lots of the layout is handled automatically by QtWidgets.QMainWindow.
        """

        self.setWindowTitle("OpenEP: The open-source solution for electrophysiology data analysis")

    def _add_menu_bar(self):
        """
        Add a menu bar to the GUI.

        The file menu has options for loading either an OpenEP dataset or an openCARP one, as well as
        loading auxillary data into an existing openCARP system.

        The view menu has options for creating a new 3D-viewer for an existing system.
        """

        self.menubar = self.menuBar()
        file_menu = self.menubar.addMenu("File")
        self.menubar.setNativeMenuBar(True)

        load_menu = QtWidgets.QMenu("Load", self)
        load_case_action = QtWidgets.QAction("OpenEP", self)
        load_case_action.triggered.connect(self._load_case)
        load_menu.addAction(load_case_action)
        load_openCARP_action = QtWidgets.QAction("openCARP", self)
        load_openCARP_action.triggered.connect(self._load_openCARP)
        load_menu.addAction(load_openCARP_action)

        # The below will be used for adding data to openCARP datasets
        self.add_data_menu = QtWidgets.QMenu("Add data to", self)

        file_menu.addMenu(load_menu)
        file_menu.addMenu(self.add_data_menu)

        # The below will be used for creating a new 3D viewer for a given system
        view_menu = self.menubar.addMenu("View")
        self.add_view_menu = QtWidgets.QMenu("Add view for", self)
        view_menu.addMenu(self.add_view_menu)
        view_menu.addSeparator()

    def _load_case(self):
        """Not yet implemented"""
        pass

    def _load_openCARP(self):
        """
        Load a set of openCARP files.

        Requires that both a *.pts and *.elem file are selected (and no other file).
        Auxillary added can later be by going to 'File > Add data to'.
        """

        dialogue = QtWidgets.QFileDialog()
        dialogue.setWindowTitle('Load a set of openCARP files')
        dialogue.setDirectory(QtCore.QDir.currentPath())
        dialogue.setFileMode(QtWidgets.QFileDialog.ExistingFiles)
        dialogue.setNameFilter("openCARP files (*.pts *.elem)")

        error = QtWidgets.QMessageBox()
        error.setIcon(QtWidgets.QMessageBox.Critical)
        error.setText("File selection error")
        error.setInformativeText(
            "Please select a single set of openCARP files (one each of *.pts, *.elem)."
        )
        error.setWindowTitle("Error")

        if dialogue.exec_():

            filenames = dialogue.selectedFiles()
            extensions = [pathlib.Path(f).suffix for f in filenames]

            if len(filenames) != 2:
                error.exec_()
            elif (".pts" not in extensions) or (".elem" not in extensions):
                error.exec_()
                return

            self._initialise_openCARP(
                points=filenames[extensions.index(".pts")],
                indices=filenames[extensions.index(".elem")],
            )

    def _initialise_openCARP(self, points, indices):
        """
        Initialise data from an openCARP simulation.

        A new system is created. If no other systems exist, this system is made the active
        system.

        Options for adding auxillary data to the system and opening a new 3D-viewer are added
        to the 'File > Add data to' and the 'View > Add view for' menus, respectively.
        """

        carp = openep.load_openCARP(
            points=points,
            indices=indices,
        )

        new_system = System(
            name=self._system_counter,
            folder=pathlib.Path(points).parent.resolve(),
            type="openCARP",
            data=carp,
        )

        self.systems.append(new_system)
        self._system_counter += 1
        self._active_system = new_system.name

        # We need to dynamically add options for loading data/creating 3D viewers to the menu
        add_data_to_system_menu = QtWidgets.QMenu(points, self)
        add_unipolar_action = QtWidgets.QAction("Unipolar electrograms", self)
        add_unipolar_action.triggered.connect(lambda: self._add_data_to_openCARP(new_system))
        add_data_to_system_menu.addAction(add_unipolar_action)
        self.add_data_menu.addMenu(add_data_to_system_menu)

        add_view_for_system_action = QtWidgets.QAction(points, self)
        add_view_for_system_action.triggered.connect(lambda: self.add_view(new_system))
        self.add_view_menu.addAction(add_view_for_system_action)

    def _add_data_to_openCARP(self, system):
        """
        Add data into an CARPData object.

        If there are no plotters for this system, a new 3D viewer will automatically be created.
        """

        dialogue = QtWidgets.QFileDialog()
        dialogue.setWindowTitle('Add an openCARP data file')
        dialogue.setDirectory(QtCore.QDir.currentPath())
        dialogue.setFileMode(QtWidgets.QFileDialog.ExistingFile)
        dialogue.setNameFilter("openCARP data file (*.dat)")

        if dialogue.exec_():

            filename = dialogue.selectedFiles()[0]
            system.data.add_data('unipolar egm', filename)
            if len(system.plotters) == 0:
                self.add_view(system)

    def add_view(self, system):
        """
        Create a new CustomDockWidget for the given system (i.e. open a new 3D viewer).

        By default, the bipolar voltage is drawn.
        """

        if system.data.electric is None:
            error = QtWidgets.QMessageBox()
            error.setIcon(QtWidgets.QMessageBox.Critical)
            error.setText("No Data Error")
            error.setInformativeText(
                "Please load data into the system before creating a new view.\n"
                "File > Add data to"
            )
            error.setWindowTitle("Error")
            error.exec_()
            return

        dock, plotter = system.create_dock()
        mesh = system.data.create_mesh()
        add_mesh_kws = system._create_default_kws()
        free_boundaries = openep.mesh.get_free_boundaries(mesh)

        index = len(system.plotters)
        plotter.lower_limit.returnPressed.connect(lambda: self.update_colourbar_limits(system, index=index))
        plotter.upper_limit.returnPressed.connect(lambda: self.update_colourbar_limits(system, index=index))

        system.docks.append(dock)
        system.plotters.append(plotter)
        system.meshes.append(mesh)
        system.add_mesh_kws.append(add_mesh_kws)
        system.free_boundaries.append(free_boundaries)

        # TODO: draw active scalars rather than bipolar_egm
        self.draw_map(
            mesh=system.meshes[index],
            plotter=system.plotters[index],
            data=system.data.electric.bipolar_egm.voltage,
            add_mesh_kws=system.add_mesh_kws[index],
        )

        self.addDockWidget(Qt.LeftDockWidgetArea, dock)
        self.add_border_to_active_plotters()

    def add_border_to_active_plotters(self):
        """Give all plotters of the active system a red border and all other plotters no border."""

        for system in self.systems:
            if system.name == self._active_system:
                for plotter in system.plotters:
                    plotter.renderer.add_border(color="red", width=5)
            else:
                for plotter in system.plotters:
                    plotter.renderer.add_border(color=None, width=5)

    def update_colourbar_limits(self, system, index):
        """Update colourbar limits for a given plotter.

        Args:
            system (System): system containing the plotter whose colourbar
                limits will be changed
            index (int): index of plotter to modify
        """

        if not system.plotters[index].lower_limit.text().strip():
            pass
        elif not system.plotters[index].upper_limit.text().strip():
            pass

        lower_limit = float(system.plotters[index].lower_limit.text())
        upper_limit = float(system.plotters[index].upper_limit.text())
        system.add_mesh_kws[index]["clim"] = [lower_limit, upper_limit]

        # TODO: draw the active scalars rather than bipolar voltage
        self.draw_map(
            system.meshes[index],
            system.plotters[index],
            system.data.electric.bipolar_egm.voltage,
            system.add_mesh_kws[index],
            system.free_boundaries[index],
        )


    def draw_map(self, mesh, plotter, data, add_mesh_kws, free_boundaries=None):
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

        if free_boundaries is not None:
            plotter = openep.draw.draw_free_boundaries(
                free_boundaries,
                colour="black",
                width=5,
                plotter=plotter,
        )

def main():

    # Create an instance of Qapplication
    app = QtWidgets.QApplication(sys.argv)
    app.setWindowIcon(QtGui.QIcon(openep.view.images.LOGO))

    # Create an instance of GUI
    window = OpenEpGUI()
    window.showMaximized()

    sys.exit(app.exec_())


@attrs(auto_attribs=True, auto_detect=True)
class System:
    """Class for storing data a single system (either an OpenEP or openCARP dataset).

    Args:
        name (str): Label for the system
        folder (str): Directory containing the file(s) for the system
        type (str): System type - either OpenEP or openCARP
        data (Union[Case, CARPData]): Data for the system, including e.g. points and triangles
            for creating a mesh.
    """
    name: str
    folder: str
    type: str
    data: Union[openep.data_structures.case.Case, openep.data_structures.openCARP.CARPData]

    def __attrs_post_init__(self):

        self.docks = []
        self.plotters = []
        self.meshes = []
        self.free_boundaries = []
        self.add_mesh_kws = []

    def __repr__(self):
        return f"{self.type} dataset {len(self.data.points)} mapping sites."

    def create_dock(self):
        """Create a new dockable pyvista-qt plotter for rendering 3D maps.

        This can be used for projecting scalar properties onto the 3D
        surface.
        """

        dock  = openep.view.custom_widgets.CustomDockWidget(f"{self.name}: Voltage")
        plotter = openep.view.plotters.create_plotter()

        plotter_layout = QtWidgets.QVBoxLayout(plotter)

        # Widgets for setting the colourbar limits
        colour_bar_layout = QtWidgets.QFormLayout()
        colour_bar_layout.setFormAlignment(QtCore.Qt.AlignLeft) 
        colour_bar_selections = QtWidgets.QHBoxLayout()

        lower_limit = QtWidgets.QLineEdit()
        lower_limit.setMinimumWidth(50)
        lower_limit.setText("0")
        lower_limit.setPlaceholderText("lower")
        lower_limit.setValidator(QtGui.QIntValidator(bottom=0.0))
        colour_bar_selections.addWidget(lower_limit)

        upper_limit = QtWidgets.QLineEdit()
        upper_limit.setMinimumWidth(50)
        upper_limit.setText("5")
        upper_limit.setPlaceholderText("upper")
        lower_limit.setValidator(QtGui.QIntValidator(bottom=0.0))
        colour_bar_selections.addWidget(upper_limit)

        # Add layouts
        colour_bar_layout.addRow("Colourbar limits:", colour_bar_selections)
        plotter_layout.addLayout(colour_bar_layout)

        # Set widget
        dock.setWidget(plotter)
        dock.setAllowedAreas(Qt.AllDockWidgetAreas)

        plotter.lower_limit = lower_limit
        plotter.upper_limit = upper_limit

        return dock, plotter

    def _create_default_kws(self):
        return {
            "clim": [0, 5],
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
