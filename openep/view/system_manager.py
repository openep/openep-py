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
Class and functions for managing OpenEP and openCARP systems loaded into the GUI.
"""
from typing import Union
from attr import attrs

from PyQt5 import QtGui, QtWidgets
from PyQt5.QtCore import Qt
import numpy as np

import openep
import openep.view.custom_widgets
import openep.view.plotters


@attrs(auto_attribs=True, auto_detect=True)
class System:
    """Class for storing data of a single system (either an OpenEP or openCARP dataset).

    Args:
        name (str): Label for the system
        basename (str): basename of the file(s) for the system
        type (str): System type - either OpenEP or openCARP
        data (Union[Case, CARPData]): Data for the system, including e.g. points and triangles
            for creating a mesh.
    """
    name: str
    basename: str
    type: str
    data: openep.data_structures.case.Case

    def __attrs_post_init__(self):

        self.docks = []
        self.plotters = []
        self.meshes = []
        self.free_boundaries = []
        self.add_mesh_kws = []
        self.scalar_fields = None

    def __repr__(self):
        return f"{self.type} dataset {len(self.data.points)} mapping sites."

    def _determine_available_fields(self):
        """Create a dictionary of scalar fields that can be used to colour the 3D map."""

        if self.type == 'OpenEP':
            self._determine_available_OpenEP_fields()
        elif self.type == 'openCARP':
            self._determine_available_openCARP_fields()

    def _determine_available_OpenEP_fields(self):
        """Create a dictionary of scalar fields that can be used to colour the 3D map."""

        self.scalar_fields = {}

        # Add scalar fields that are interpolated by openep
        if self.data.electric.bipolar_egm.egm is not None:
            field = openep.case.interpolate_voltage_onto_surface(
                self.data,
                max_distance=None,
                buffer=0,
            )
            self.scalar_fields['Bipolar voltage'] = field

        if self.data.electric.unipolar_egm.egm is not None:
            field = openep.case.interpolate_voltage_onto_surface(
                self.data,
                max_distance=None,
                buffer=0,
                bipolar=False,
            )
            self.scalar_fields['Unipolar voltage'] = field

        # Add scalar fields that are calculate by the clinical mapping system
        clinical_field_names = [
            'Clinical bipolar voltage',
            'Clinical unipolar voltage',
            'Clinical LAT',
            'Clinical force',
            'Clinical impedance ',
        ]
        clinical_fields = [
            self.data.fields.bipolar_voltage,
            self.data.fields.unipolar_voltage,
            self.data.fields.local_activation_time,
            self.data.fields.force,
            self.data.fields.impedance,
        ]

        for field_name, field in zip(clinical_field_names, clinical_fields):
            
            if np.isnan(field).all():
                continue

            self.scalar_fields[field_name] = field

    def _determine_available_openCARP_fields(self):
        """Create a dictionary of scalar fields that can be used to colour the 3D map."""

        self.scalar_fields = {}

        if self.data.electric.bipolar_egm.voltage is not None:
            self.scalar_fields['Bipolar voltage'] = self.data.electric.bipolar_egm.voltage

        if self.data.electric.unipolar_egm.voltage is not None:
            self.scalar_fields['Unipolar voltage'] = self.data.electric.unipolar_egm.voltage[:, 0]

    def create_dock(self):
        """Create a new dockable pyvista-qt plotter for rendering 3D maps.

        This can be used for projecting scalar properties onto the 3D
        surface.
        """

        if self.data is None or self.data.electric is None:
            return

        if self.scalar_fields is None:
            self._determine_available_fields()

        dock = openep.view.custom_widgets.CustomDockWidget("Temporary title")  # this will be changed
        dock.main = QtWidgets.QMainWindow()
        plotter = openep.view.plotters.create_plotter()
        mesh = self.data.create_mesh()

        self._add_scalars_to_mesh(mesh)

        # TODO: We're currently using line edits to set colourbar limits - look into RangeSlider
        plotter_layout = QtWidgets.QVBoxLayout(plotter)
        plotter.colour_bar_layout, plotter.lower_limit, plotter.upper_limit = self._create_colourbar_layout()
        plotter.opacity_layout, plotter.opacity = self._create_opacity_layout()
        control_layout = self._create_control_layout(plotter=plotter)
        plotter_layout.addLayout(control_layout)
        plotter_layout.addStretch()

        # The dock is set to have bold font (so the title stands out)
        # But all other widgets should have normal weighted font
        main_font = QtGui.QFont()
        main_font.setBold(False)
        dock.main.setFont(main_font)

        # Create a menubar for this dock
        # From here we can control e.g. the scalar values projected onto the surface
        dock.main.menubar = dock.main.menuBar()
        dock.main.menubar.setNativeMenuBar(False)
        dock.main.menubar

        self._add_field_menu(dock, plotter)

        # Set widget
        dock.main.setCentralWidget(plotter)
        dock.setWidget(dock.main)
        dock.setAllowedAreas(Qt.AllDockWidgetAreas)

        return dock, plotter, mesh

    def _add_scalars_to_mesh(self, mesh):
        """Add all scalar fields to the mesh point data"""

        if self.scalar_fields is None:
            return

        for field_name, field in self.scalar_fields.items():
            mesh.point_data[field_name] = field

    def _create_colourbar_layout(self):
        """Create a layout with widgets for setting the limits of the colourbar."""

        colour_bar_layout = QtWidgets.QHBoxLayout()

        colour_bar_text = QtWidgets.QLabel("Colourbar limits:")
        colour_bar_text.setMinimumWidth(100)
        colour_bar_text.setMaximumWidth(200)
        colour_bar_text.setStyleSheet('QLabel {border: 0px; background-color: #d8dcd6;}')
        colour_bar_layout.addWidget(colour_bar_text)

        lower_limit = QtWidgets.QLineEdit()
        lower_limit.setFixedWidth(40)
        lower_limit.setText("0")
        lower_limit.setPlaceholderText("lower")
        lower_limit.setValidator(QtGui.QDoubleValidator(bottom=0.0))
        colour_bar_layout.addWidget(lower_limit)

        upper_limit = QtWidgets.QLineEdit()
        upper_limit.setFixedWidth(40)
        upper_limit.setText("5")
        upper_limit.setPlaceholderText("upper")
        lower_limit.setValidator(QtGui.QDoubleValidator(bottom=0.01))
        colour_bar_layout.addWidget(upper_limit)

        colour_bar_layout.addStretch()

        return colour_bar_layout, lower_limit, upper_limit

    def _create_opacity_layout(self):
        """Create a layout with widgets for setting the opacity of the mesh."""

        opacity_layout = QtWidgets.QHBoxLayout()
        opacity_layout.addStretch

        opacity_text = QtWidgets.QLabel("Opacity:")
        opacity_text.setMinimumWidth(50)
        opacity_text.setMaximumWidth(120)
        opacity_text.setStyleSheet('QLabel {border: 0px; background-color: #d8dcd6;}')
        opacity_layout.addWidget(opacity_text)

        opacity_selector = QtWidgets.QDoubleSpinBox()
        opacity_selector.setFixedWidth(60)
        opacity_selector.setRange(0.0, 1.0)
        opacity_selector.setSingleStep(0.05)
        opacity_selector.setDecimals(2)
        opacity_selector.setWrapping(False)
        opacity_selector.setValue(1)
        opacity_selector.textFromValue(1)
        opacity_layout.addWidget(opacity_selector)

        return opacity_layout, opacity_selector

    def _create_control_layout(self, plotter):
        """Add the colourbar and opacity layout to a single horizontal layout"""

        control_layout = QtWidgets.QHBoxLayout()
        control_layout.addLayout(plotter.colour_bar_layout)
        control_layout.addLayout(plotter.opacity_layout)

        return control_layout

    def _add_field_menu(self, dock, plotter):
        """Add a Field menu to the menubar. This is used for selecting the scalar field to project onto the surface."""

        field_menu = dock.main.menubar.addMenu("Field")
        field_group = QtWidgets.QActionGroup(dock.main)

        plotter.scalar_field_actions = {}
        for field_name in self.scalar_fields:

            dock.setWindowTitle(f"{self.name}: {field_name}")
            action = QtWidgets.QAction(field_name, dock.main, checkable=True)
            action.setChecked(False)

            field_menu.addAction(action)
            field_group.addAction(action)
            plotter.scalar_field_actions[field_name] = action

    def _create_default_kws(self):
        """Create a dictionary of keyword arguments to be passed to openep.draw.draw_map"""

        return {
            "clim": [0, 5],
            "scalar_bar_args": {
                "title": "Voltage (mV)",
            }
        }
