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

        self.setDockOptions(self.GroupedDragging | self.AllowTabbedDocks | self.AllowNestedDocks)
        self.setTabPosition(Qt.AllDockWidgetAreas, QtWidgets.QTabWidget.North)

    def _init_data(self):
        """
        Set default values for variables that can later be set by the user.
        """

        # Initially, display first electrogram only
        self.egm_points = np.array([0], dtype=int)

        # initial bipolar voltage colourbar limits (for when new windows are opened)
        self._initial_colourbar_limits = (0, 2)

        # We need a list of all cases currently loaded
        # Once a case is loaded, we should keep track of:
        # * it's name (defaults to monotonically increasing integers)
        # * it's folder
        # * it's type (OpenEP, openCARP)
        # * it's data - either a Case or CARPData object
        # * a list of it's dock widgets (one plotter per dock)
        # i.e. {'name': str, 'folder': str, 'type': str, 'data': Union[Case, CARPData], 'docks': list[DockWidget]}
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
        """Add a menu bar to the GUI."""

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

        file_menu.addMenu(load_menu)

    def _load_case(self):
        pass

    def _load_openCARP(self):
        """Load a set of openCARP files."""

        dialogue = QtWidgets.QFileDialog()
        dialogue.setWindowTitle('Load a set of openCARP files')
        dialogue.setDirectory(QtCore.QDir.currentPath())
        dialogue.setFileMode(QtWidgets.QFileDialog.ExistingFiles)
        dialogue.setNameFilter("openCARP files (*.pts *.elem *.dat)")

        error = QtWidgets.QMessageBox()
        error.setIcon(QtWidgets.QMessageBox.Critical)
        error.setText("File selection error")
        error.setInformativeText(
            "Please select a single set of openCARP files (one each of *.pts, *.elem, *.dat)."
        )
        error.setWindowTitle("Error")

        if dialogue.exec_():

            filenames = dialogue.selectedFiles()
            extensions = [pathlib.Path(f).suffix for f in filenames]

            if len(filenames) != 3:
                error.exec_()
            elif (".pts" not in extensions) or (".elem" not in extensions) or (".dat" not in extensions):
                error.exec_()

            self._initialise_openCARP(
                points=filenames[extensions.index(".pts")],
                indices=filenames[extensions.index(".elem")],
                data=filenames[extensions.index(".dat")],
            )

    def _initialise_openCARP(self, points, indices, data):
        """
        Initialise data from an openCARP simulation.
        """

        carp = openep.load_openCARP(
            points=points,
            indices=indices,
            unipolar_egm=data,
        )

        new_system = System(
            name=self._system_counter,
            folder=pathlib.Path(points).parent.resolve(),
            type="openCARP",
            data=carp,
        )

        # TODO: Creation of new dock etc. for a given system should be put into a separate function.
        #       This will allow us to create new docks etc. for a system that already exists.
        #       e.g. add a `Add View` menu, can select the system from the list for which we want a
        #       new view (i.e. a new dock widgets and associated plotter etc.)
        dock, plotter = new_system.create_dock()
        mesh = carp.create_mesh()
        add_mesh_kws = new_system._create_default_kws()
        free_boundaries = openep.mesh.get_free_boundaries(mesh)

        # This is the first plotter added to the system, so the index is 0
        plotter.lower_limit.returnPressed.connect(lambda: self.update_colourbar_limits(new_system, index=0))
        plotter.upper_limit.returnPressed.connect(lambda: self.update_colourbar_limits(new_system, index=0))

        new_system.docks.append(dock)
        new_system.plotters.append(plotter)
        new_system.meshes.append(mesh)
        new_system.add_mesh_kws.append(add_mesh_kws)
        new_system.free_boundaries.append(free_boundaries)

        # TODO: draw active scalars rather than bipolar_egm
        self.draw_map(
            mesh=new_system.meshes[0],
            plotter=new_system.plotters[0],
            data=new_system.data.bipolar_egm,
            add_mesh_kws=new_system.add_mesh_kws[0],
        )

        self.addDockWidget(Qt.LeftDockWidgetArea, dock)
        self.systems.append(new_system)
        self._system_counter += 1
        self._active_system = new_system.name
        self.add_border_to_active_plotters()

    def add_border_to_active_plotters(self):

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
        print(system, system.add_mesh_kws[index]['clim'])
        self.draw_map(
            system.meshes[index],
            system.plotters[index],
            system.data.bipolar_egm,
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
        return f"{self.type} dataset {len(self.data.unipolar_egm)} mapping sites."

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

# Will we always have unipolar voltages form openCARP?
# We need a way to optionally take, unipolar, bipolar, local activation time,
# conduction velocity, etc.
# These can be taken as optional arguments to CARPData. Do not automatically generate the bipolar data.
# Instead, allow the user to call CARPData.bipolar_from_unipolar()
# Then we need a way to calculate e.g. local activation times, conduction velocity

# We also need a way to load data into openCARP systems from the GUI. Perhaps from the menu
# File > Add data > (select system) > (select data type, e.g. unipolar) > open file dialogue

# How can we allow other simulation formats? e.g. CHASTE

