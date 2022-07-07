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
        bipolar_voltage (np.ndarray): array of shape N
        unipolar_voltage (np.ndarray): array of shape N
        local_activation_time (np.ndarray): array of shape N
        impedance (np.ndarray): array of shape N
        force (np.ndarray): array of shape N
    """

    bipolar_voltage: np.ndarray
    unipolar_voltage: np.ndarray
    local_activation_time: np.ndarray
    impedance: np.ndarray
    force: np.ndarray
    thickness: np.ndarray = None

    def __attrs_post_init__(self):
        if self.bipolar_voltage is not None and self.thickness is None:
            self.thickness = np.full_like(self.bipolar_voltage, fill_value=np.NaN, dtype=float)

    def __repr__(self):
        return f"fields: {tuple(self.__dict__.keys())}"

    def __getitem__(self, field):
        try:
            return self.__dict__[field]
        except KeyError as e:
            raise ValueError(f"There is no field '{field}'") from e

    def __contains__(self, field):
        return field in self.__dict__.keys()


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
        return np.array([]), np.array([]), empty_fields(size=0)

    points = surface_data['triRep']['X'].astype(float)
    indices = surface_data['triRep']['Triangulation'].astype(int)

    if surface_data['act_bip'].size == 0:
        size = points.size // 3
        fields = empty_fields(size=size)
        return points, indices, fields

    local_activation_time, bipolar_voltage = surface_data['act_bip'].T.astype(float)
    unipolar_voltage, impedance, force = surface_data['uni_imp_frc'].T.astype(float)

    try:
        thickness = surface_data['thickness'].astype(float)
    except KeyError as e:
        thickness = np.full_like(local_activation_time, fill_value=np.NaN, dtype=float)

    fields = Fields(
        bipolar_voltage,
        unipolar_voltage,
        local_activation_time,
        impedance,
        force,
        thickness
    )

    return points, indices, fields


def empty_fields(size=0):
    """Create an empty Fields object with empty numpy arrays.

    Returns:
        fields (Field): Field object that contains numpy arrays of the various
            scalar fields
    """

    local_activation_time = np.full(size, fill_value=np.NaN, dtype=float)
    bipolar_voltage = np.full(size, fill_value=np.NaN, dtype=float)
    unipolar_voltage = np.full(size, fill_value=np.NaN, dtype=float)
    impedance = np.full(size, fill_value=np.NaN, dtype=float)
    force = np.full(size, fill_value=np.NaN, dtype=float)
    thickness = np.full(size, fill_value=np.NaN, dtype=float)

    fields = Fields(
        bipolar_voltage,
        unipolar_voltage,
        local_activation_time,
        impedance,
        force,
        thickness
    )

    return fields
