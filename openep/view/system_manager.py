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

    def create_mesh(self):
        """Create a new mesh from the system data"""

        mesh = self.data.create_mesh()
        return self._add_scalars_to_mesh(mesh)

    def _add_scalars_to_mesh(self, mesh):
        """Add all scalar fields to the mesh point data"""

        if self.scalar_fields is None:
            return

        for field_name, field in self.scalar_fields.items():
            mesh.point_data[field_name] = field

        return mesh

    def _create_default_kws(self):
        """Create a dictionary of keyword arguments to be passed to openep.draw.draw_map"""

        return {
            "clim": [0, 5],
            "scalar_bar_args": {
                "title": "Voltage (mV)",
            }
        }
