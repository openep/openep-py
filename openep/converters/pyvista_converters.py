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
PyVista converters --- :mod:`openep.converters.pyvista_converters`
==================================================================

This module contains functions to convert OpenEP datasets to/from PyVista PolyData objects.

"""

import numpy as np
import pyvista

from ..data_structures.surface import Fields
from ..data_structures.electric import Electric
from ..data_structures.ablation import Ablation
from ..data_structures.case import Case

__all__ = ["to_pyvista", "from_pyvista"]


def from_pyvista(
    mesh: pyvista.PolyData,
    name: str = None,
    scale_points: float = 1,
) -> Case:
    """Convert a pyvista mesh to an OpenEP dataset.

    Args:
        mesh (pyvista.PolyData): Mesh to convert.
        name (str, optional): Name of the dataset. If None, will be set to the memory address of the case object.
        scale_points (float, optional): Scale the point positions by this number. Useful to e.g. scale
            the units to be in mm rather than micrometre. Default is 1.

    Returns:
        case (openep.Case): OpenEP case object.

    Note
    ----
    All other attributes of the Case object will be empty numpy arrays.

    """

    case = Case(
        name=name,
        points=mesh.points * scale_points,
        indices=np.array(mesh.faces).reshape(mesh.n_cells, 4)[:, 1:],
        fields=Fields.from_pyvista(mesh),
        electric = Electric(),
        ablation = Ablation(),
        notes = np.asarray([], dtype=object),
    )

    return case


def to_pyvista(
    case: Case,
    add_field_data: bool = False,
) -> pyvista.PolyData:
    """Convert an OpenEP dataset to a PyVista.Polydata mesh.

    Args:
        case (openep.Case): Case to convert to a PyVista mesh.
        add_field_data (bool, optional). If True, add the data from case.field
            as point data in the new mesh (mesh.point_data).

    Returns:
        mesh (pyvista.PolyData): 
    """

    mesh = case.create_mesh()

    if add_field_data:
        mesh.point_data["Bipolar voltage"] = case.fields.bipolar_voltage
        mesh.point_data["Unipolar voltage"] = case.fields.unipolar_voltage
        mesh.point_data["LAT"] = case.fields.local_activation_time
        mesh.point_data["Impedance"] = case.fields.impedance
        mesh.point_data["Force"] = case.fields.force
        mesh.point_data["Conduction velocity"] = case.fields.conduction_velocity
        mesh.point_data["Divergence"] = case.fields.cv_divergence

    return mesh
