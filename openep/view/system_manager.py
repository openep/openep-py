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

import numpy as np
import pyvista

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
        case (Case): Data for the system, including e.g. points and triangles
            for creating a mesh.
    """
    name: str
    basename: str
    type: str
    case: openep.data_structures.case.Case

    def __attrs_post_init__(self):

        self.docks = []
        self.plotters = []
        self.meshes = []
        self.mapping_points_meshes = []
        self.surface_projected_discs_meshes = []
        self.free_boundaries = []
        self.add_mesh_kws = []
        self.scalar_fields = None
        self.egm_points = np.asarray([0], dtype=int)  # default selected electrograms

    def __repr__(self):
        return f"{self.type} dataset {len(self.case.points)} mapping sites."

    def _determine_available_fields(self):
        """Create a dictionary of scalar fields that can be used to colour the 3D map."""

        # TODO: Don't distinguish between Clinical and non-clinical data. Use clinical by default. If
        #       values are interpolated from mapping points onto the surface, then overwrite the
        #       default clinical values. Add an option to reset the interpolated values to the original
        #       clinical ones. This can be done by simply reading the file.

        if self.type == 'OpenEP':
            self._determine_available_OpenEP_fields()
        elif self.type == 'openCARP':
            self._determine_available_openCARP_fields()
        else:
            raise ValueError("self.type must be one of: OpenEP; openCARP")

    def _determine_available_OpenEP_fields(self):
        """Create a dictionary of scalar fields that can be used to colour the 3D map."""

        self.scalar_fields = {}

        # Add scalar fields that are interpolated by openep
        if self.case.electric.bipolar_egm.egm is not None:
            field = openep.case.interpolate_voltage_onto_surface(
                self.case,
                max_distance=None,
                buffer=0,
            )
            self.scalar_fields['Bipolar voltage'] = field

        if self.case.electric.unipolar_egm.egm is not None:
            field = openep.case.interpolate_voltage_onto_surface(
                self.case,
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
            self.case.fields.bipolar_voltage,
            self.case.fields.unipolar_voltage,
            self.case.fields.local_activation_time,
            self.case.fields.force,
            self.case.fields.impedance,
        ]

        for field_name, field in zip(clinical_field_names, clinical_fields):

            if np.isnan(field).all():
                continue

            self.scalar_fields[field_name] = field

    def _determine_available_openCARP_fields(self):
        """Create a dictionary of scalar fields that can be used to colour the 3D map."""

        self.scalar_fields = {}

        if self.case.electric.bipolar_egm.voltage is not None:
            self.scalar_fields['Bipolar voltage'] = self.case.electric.bipolar_egm.voltage

        if self.case.electric.unipolar_egm.voltage is not None:
            self.scalar_fields['Unipolar voltage'] = self.case.electric.unipolar_egm.voltage[:, 0]

    def create_mesh(self):
        """Create a new mesh from the system data"""

        mesh = self.case.create_mesh()
        return self._add_scalars_to_mesh(mesh)

    def _add_scalars_to_mesh(self, mesh):
        """Add all scalar fields to the mesh point data"""

        if self.scalar_fields is None:
            return

        for field_name, field in self.scalar_fields.items():
            mesh.point_data[field_name] = field

        return mesh

    def create_selected_point_mesh(self, point):
        """Create a mesh that will be used for highlighting the selected point"""

        point_mesh = pyvista.PolyData(point)
        point_geometry = pyvista.Sphere(theta_resolution=20, phi_resolution=20)

        factor = 2 if self.type == "OpenEP" else 2000
        glyphed_mesh = point_mesh.glyph(scale=False, factor=factor, geom=point_geometry)

        return glyphed_mesh

    def create_mapping_points_mesh(self):
        """Create a mesh that contains the mapping points."""

        mapping_points_centered = self.case.electric.bipolar_egm.points - self.case._mesh_center
        mapping_points_mesh = pyvista.PolyData(mapping_points_centered)
        point_geometry = pyvista.Sphere(theta_resolution=8, phi_resolution=8)
        
        factor = 1.5 if self.type == "OpenEP" else 1200
        glyphed_mesh = mapping_points_mesh.glyph(scale=False, factor=factor, geom=point_geometry)

        return glyphed_mesh

    def create_surface_discs_mesh(self):
        """Create a mesh that contains the surface-projected mapping points"""

        mapping_points_centered = self.case.electric.bipolar_egm.points - self.case._mesh_center
        surface_points_centered = self.case.points - self.case._mesh_center

        distance_to_surface = openep.case.case_routines.calculate_distance(
            origin=mapping_points_centered,
            destination=surface_points_centered
        )
        nearest_point_indices = np.argmin(distance_to_surface, axis=1)

        projected_mapping_points = surface_points_centered[nearest_point_indices]
        projected_mapping_points_mesh = pyvista.PolyData(projected_mapping_points)

        # TODO: we're currently plotting these points as spheres. They should be discs instead, oriented
        #       along the surface of the mesh.
        disc_geometry = pyvista.Sphere(theta_resolution=8, phi_resolution=8)
        factor = 1.5 if self.type == "OpenEP" else 1200
        glyphed_mesh = projected_mapping_points_mesh.glyph(
            scale=False,
            factor=factor,
            geom=disc_geometry,
        )

        return glyphed_mesh

    def _create_default_kws(self):
        """Create a dictionary of keyword arguments for plotting meshes, points, and surface-projected discs."""

        add_mesh_kws = {
            "pickable": False,  # never set to True
            "clim": [0, 5],
            "name": "Surface",
            "opacity": 1,
            "scalar_bar_args": {
                "title": "Voltage (mV)",
                "color": "white",
            }
        }

        add_points_kws = {
            "pickable": False,  # only set to True for the primary viewer
            "name": "Mapping points",
            "opacity": 1,
            "show_scalar_bar": False,
            "smooth_shading": True,
            "lighting": True,
        }

        add_discs_kws = {
            "pickable": False,  # only set to True for the primary viewer
            "name": "Surface-projected mapping points",
            "opacity": 1,
            "show_scalar_bar": False,
            "smooth_shading": True,
            "lighting": True,
        }

        return add_mesh_kws, add_points_kws, add_discs_kws


class SystemManager:
    """Keep track of all systems loaded into the GUI"""

    def __init__(self):

        # We need to keep track of all cases currently loaded (both OpenEP and openCARP)
        # Each system will keep track of its own docks/plotters/meshes
        self.systems = {}
        self.system_counter = 0  # we need a unique ID for each new system (IDs of deleted systems should not be reused)

        # We will use this to determine which egms to plot
        # All windows belonging to the main case have their menubars coloured blue
        self.active_system = None

    def add_system(self, basename, type, case):
        """Add a system to the dictionary of systems.

        If this is the first system, set it to be the active system.

        Args:
            basename (str): basename of the file(s) for the system
            type (str): System type - either OpenEP or openCARP
            case (Case): Data for the system, including e.g. points and triangles
                for creating a mesh.

        Returns:
            system (System): newly created system
        """

        # we need a unique key for each system
        while self.system_counter in self.systems:
            self.system_counter += 1

        name = str(self.system_counter)
        system = System(
            name=name,
            basename=basename,
            type=type,
            case=case,
        )

        self.systems[name] = system
        self.system_counter += 1
        self.active_system = self.active_system or system

        system._determine_available_fields()

        return system

    def update_system_name(self, system, new_name):
        """Change the name of a system. It must be unique - two systems cannot have the same name.py

        The name is updated in-place.

        Args:
            system (System): system whose name will be changed
            new_name (str): new name for the system
        """

        if new_name in self.systems:
            raise KeyError("System names must be unique")

        self.systems[new_name] = self.systems.pop(system.name)
        self.systems[new_name].name = new_name

    def check_available_electrograms(self):
        """Check which electrogram type(s) are present in the active system.

        Returns:
            has_reference (bool): If True, the active system has reference electrograms
            has_bipolar (bool): If True, the active system has bipolar electrograms
            has_unipolar (bool): If True, the active system has unipolar electrograms
        """

        electric = self.active_system.case.electric

        has_reference = True if electric.reference_egm.egm is not None else False
        has_bipolar = True if electric.bipolar_egm.egm is not None else False
        has_unipolar = True if electric.unipolar_egm.egm is not None else False

        return has_reference, has_bipolar, has_unipolar

    def electrogram_times(self):
        """Generate the electrogram times based on the number of values in each electrogram.

        Warning
        -------
        This assumes that the number of values in each type of electrogram (reference,
        bipolar, unipolar) is the same.
        """

        electric = self.active_system.case.electric
        egm_times = None

        try:
            return np.arange(electric.reference_egm.egm.shape[1])
        except AttributeError:
            pass

        try:
            return np.arange(electric.biipolar_egm.egm.shape[1])
        except AttributeError:
            pass

        try:
            return np.arange(electric.unipolar_egm.egm.shape[1])
        except AttributeError:
            pass

        return egm_times
