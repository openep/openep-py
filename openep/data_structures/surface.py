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

"""Module containing classes for storing surface data of a mesh."""

from attr import attrs
import numpy as np

__all__ = []


@attrs(auto_attribs=True, auto_detect=True)
class Fields:
    """
    Class for storing information about the surface of a mesh

    Args:
        bipolar_voltage (np.ndarray): array of shape N_points
        unipolar_voltage (np.ndarray): array of shape N_points
        local_activation_time (np.ndarray): array of shape N_points
        impedance (np.ndarray): array of shape N_points
        force (np.ndarray): array of shape N_points
        region (np.ndarray): array of shape N_cells
        longitudinal_fibres (np.ndarray): array of shape N_cells x 3
        transverse_fibres (np.ndarray): array of shape N_cells x 3
        pacing_site (np.ndarray): array of shape N_points
    """

    bipolar_voltage: np.ndarray = None
    unipolar_voltage: np.ndarray = None
    local_activation_time: np.ndarray = None
    impedance: np.ndarray = None
    force: np.ndarray = None
    thickness: np.ndarray = None
    cell_region: np.ndarray = None
    longitudinal_fibres: np.ndarray = None
    transverse_fibres: np.ndarray = None
    pacing_site: np.ndarray = None
    conduction_velocity: np.ndarray = None
    cv_divergence: np.ndarray = None
    mesh_normals : np.ndarray = None

    def __repr__(self):
        return f"fields: {tuple(self.__dict__.keys())}"

    def __getitem__(self, field):
        try:
            return self.__dict__[field]
        except KeyError as e:
            raise ValueError(f"There is no field '{field}'.") from e

    def __setitem__(self, field, value):
        if field not in self.__dict__.keys():
            raise ValueError(f"'{field}' is not a valid field name.")
        self.__dict__[field] = value

    def __iter__(self):
        return iter(self.__dict__.keys())

    def __contains__(self, field):
        return field in self.__dict__.keys()

    def copy(self):
        """Create a deep copy of Fields"""

        fields = Fields()
        for field in self:
            if self[field] is None:
                continue
            fields[field] = np.array(self[field])

        return fields

    @classmethod
    def from_pyvista(cls, mesh):
        """Create a Fields object from a pyvista.PolyData mesh.

        Args:
            mesh (pyvista.PolyData): mesh from which Fields will be created.
        """

        fields = cls()
        for point_data in mesh.point_data:
            if point_data not in fields:
                continue
            fields[point_data] = np.asarray(mesh.point_data[point_data])
        for cell_data in mesh.cell_data:
            if cell_data not in fields:
                continue
            fields[cell_data] = np.asarray(mesh.cell_data[cell_data])

        return fields


def extract_surface_data(surface_data):
    """Extract surface data from a dictionary.

    Args:
        surface_data (dict): Dictionary containing numpy arrays that describe the
            surface of a mesh as well as scalar values (fields)

    Returns:
        points (ndarray): 3D coordinates of points on the mesh
        indices (ndarray): Indices of points that make up each face of the mesh
        fields (Field): Field object that contains numpy arrays of the various
            scalar fields
    """

    if surface_data['triRep']['X'].size == 0:
        return None, None, Fields()

    points = surface_data['triRep']['X'].astype(float)
    indices = surface_data['triRep']['Triangulation'].astype(int)

    if surface_data['act_bip'].size == 0:
        local_activation_time = None
        bipolar_voltage = None
    else:
        local_activation_time, bipolar_voltage = surface_data['act_bip'].T.astype(float)

    if isinstance(local_activation_time, np.ndarray) and all(np.isnan(local_activation_time)):
        local_activation_time = None
    if isinstance(bipolar_voltage, np.ndarray) and all(np.isnan(bipolar_voltage)):
        bipolar_voltage = None

    if surface_data['uni_imp_frc'].size == 0:
        unipolar_voltage = None
        impedance = None
        force = None
    else:
        unipolar_voltage, impedance, force = surface_data['uni_imp_frc'].T.astype(float)
        if all(np.isnan(unipolar_voltage)):
            unipolar_voltage = None
        if all(np.isnan(impedance)):
            impedance = None
        if all(np.isnan(force)):
            force = None

    try:
        thickness = surface_data['thickness'].astype(float)
    except KeyError as e:
        thickness = None

    if isinstance(thickness, np.ndarray) and thickness.size == 0:
        thickness = None

    # This is defined on a per-cell bases
    try:
        cell_region = surface_data['cell_region'].astype(int)
    except KeyError as e:
        cell_region = None

    if isinstance(cell_region, np.ndarray) and cell_region.size == 0:
        cell_region = None

    # Fibre orientation are vectors defined on a per-cell basis
    try:
        longitudinal_fibres = surface_data['fibres']['longitudinal'].astype(float)
    except:
        longitudinal_fibres = None

    if isinstance(longitudinal_fibres, np.ndarray) and longitudinal_fibres.size == 0:
        longitudinal_fibres = None

    try:
        transverse_fibres = surface_data['fibres']['transverse'].astype(float)
    except:
        transverse_fibres = None

    if isinstance(transverse_fibres, np.ndarray) and transverse_fibres.size == 0:
        transverse_fibres = None

    # Pacing site point ids (-1 is not pacing site)
    try:
        pacing_site = surface_data['pacing_site'].astype(int)
    except KeyError as e:
        pacing_site = None

    if isinstance(pacing_site, np.ndarray) and pacing_site.size == 0:
        pacing_site = None
    
    try:
        conduction_velocity = surface_data['signalMaps']['conduction_velocity_field']['value'].astype(float)
    except KeyError as e:
        conduction_velocity = None

    try:
        cv_divergence = surface_data['signalMaps']['divergence_field']['value'].astype(float)
    except KeyError as e:
        cv_divergence = None

    try:
        mesh_normals = surface_data['mesh_normals'].astype(float)
    except KeyError as e:
        mesh_normals = None

    fields = Fields(
        bipolar_voltage=bipolar_voltage,
        unipolar_voltage=unipolar_voltage,
        local_activation_time=local_activation_time,
        impedance=impedance,
        force=force,
        thickness=thickness,
        cell_region=cell_region,
        longitudinal_fibres=longitudinal_fibres,
        transverse_fibres=transverse_fibres,
        pacing_site=pacing_site,
        conduction_velocity=conduction_velocity,
        cv_divergence=cv_divergence,
        mesh_normals=mesh_normals
    )

    return points, indices, fields


def empty_fields(n_points=0, n_cells=0):
    """Create an empty Fields object with empty numpy arrays.

    Returns:
        fields (Field): Field object that contains numpy arrays of the various
            scalar fields
    """

    local_activation_time = np.full(n_points, fill_value=np.NaN, dtype=float)
    bipolar_voltage = np.full(n_points, fill_value=np.NaN, dtype=float)
    unipolar_voltage = np.full(n_points, fill_value=np.NaN, dtype=float)
    impedance = np.full(n_points, fill_value=np.NaN, dtype=float)
    force = np.full(n_points, fill_value=np.NaN, dtype=float)
    thickness = np.full(n_points, fill_value=np.NaN, dtype=float)
    cell_region = np.full(n_cells, fill_value=0, dtype=int)
    longitudinal_fibres = np.full((n_cells, 3), fill_value=np.NaN)
    transverse_fibres = np.full((n_cells, 3), fill_value=np.NaN)
    pacing_site = np.full(n_points, fill_value=-1, dtype=int)

    fields = Fields(
        bipolar_voltage,
        unipolar_voltage,
        local_activation_time,
        impedance,
        force,
        thickness,
        cell_region,
        longitudinal_fibres,
        transverse_fibres,
        pacing_site,
    )

    return fields
