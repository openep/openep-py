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
`Case` - the fundamental data structure of `openep-py`
======================================================

A `Case` stores all the information obtained during a clinical mapping
procedure.

Warning
-------
    This class should not be instatiated directly. Instead, a `Case` can be
    created using :func:`openep.io.readers.load_openep_mat`.

Tip
---
    Once you have a `Case`, you can perform analyses using the functions in
    :mod:`openep.case.case_routines`. You can also create a 3D mesh using the
    :func:`create_mesh` method.

.. autoclass:: Case
    :members: create_mesh, get_surface_data, get_field

Note
----

`Case` contains information on the 3D surface and associated scalar fields,
the electrograms, and the ablation sites. These data are, respectively,
stored as the following data structures: :class:`openep.data_structures.surface.Fields`,
:class:`openep.data_structures.electric.Electric`, and
:class:`openep.data_structures.ablation.Ablation`.


Scalar field data
-----------------

.. autoclass:: openep.data_structures.surface.Fields


Electrical data
---------------

.. autoclass:: openep.data_structures.electric.Electric

.. autoclass:: openep.data_structures.electric.Electrogram

.. autoclass:: openep.data_structures.electric.Impedance

.. autoclass:: openep.data_structures.electric.ElectricSurface

.. autoclass:: openep.data_structures.electric.Annotations


Ablation data
-------------

.. autoclass:: openep.data_structures.ablation.Ablation

.. autoclass:: openep.data_structures.ablation.AblationForce
"""

from attr import attrs
from typing import Optional, Tuple, List

import numpy as np
import pyvista

from .surface import Fields
from .electric import Electric, Electrogram, Annotations
from .ablation import Ablation
from ..case.case_routines import bipolar_from_unipolar_surface_points

__all__ = []


@attrs(auto_attribs=True, auto_detect=True)
class Case:
    """
    The fundamental OpenEP object.

    The class contains all the information on a single case exported from a clinical
    mapping system.

    Args:
        name (str): Name to assign to the dataset
        points (ndarray): 3D coordinates of points on the mesh
        indices (ndarray): Indices of points that make up each face of the mesh
        fields (Fields): Numpy arrays of the scalar fields associated with each point on
            the surface of a mesh
        electtic (Electric): Electrical data obtained during a clinical mapping procedure.
        ablation (Ablation, optional): Ablation data obtained during a clinical mapping procedure.
        notes (list, optional): Notes associated with the dataset.

    """
    name: str
    points: np.ndarray
    indices: np.ndarray
    fields: Fields
    electric: Electric
    ablation: Optional[Ablation] = None
    notes: Optional[List] = None

    def __attrs_post_init__(self):

        tmp_mesh = self.create_mesh(recenter=False)
        self._mesh_center = np.asarray(tmp_mesh.center)

    def __repr__(self):
        return f"{self.name}( nodes: {self.points.shape} indices: {self.indices.shape} {self.fields} )"

    def create_mesh(
        self,
        recenter: bool = True,
        back_faces: bool = False,
    ):
        """
        Create a new mesh object from the stored nodes and indices

        Args:
            recenter: if True, recenter the mesh to the origin
            back_faces: if True, calculate back face triangles

        Returns:
            mesh (pyvista.Polydata): a mesh created from the case's points and indices
        """

        indices = self.indices
        num_points_per_face = np.full(shape=(len(indices)), fill_value=3, dtype=int)  # all faces have three vertices
        faces = np.concatenate([num_points_per_face[:, np.newaxis], indices], axis=1)

        if back_faces:
            indices_inverted = indices[:, [0, 2, 1]]
            faces_inverted = np.concatenate([num_points_per_face[:, np.newaxis], indices_inverted], axis=1)
            faces = np.vstack(
                [faces, faces_inverted]
            )  # include each face twice for both surfaces

        mesh = pyvista.PolyData(self.points, faces.ravel())

        if recenter:
            mesh.translate(
                -np.asarray(mesh.center)
            )  # recenter mesh to origin, helps with lighting in default scene

        return mesh

    def get_surface_data(self, copy: bool = False) -> Tuple[np.ndarray, np.ndarray]:
        """
        Extract the surface data for the case.

        Args:
            copy (bool, optional): If True, a copy of the data will be returned. The default
                is False, in which case a view of the data is returned.

        Returns:
           points (ndarray): 3D coordinates of points on the mesh
           indices (ndarray): Indices of points that make up each face of the mesh
        """
        points = self.points
        indices = self.indices

        if copy:
            points = points.copy()
            indices = indices.copy()

        return points, indices

    def get_field(self, fieldname: str, copy: bool = False) -> np.ndarray:
        """
        Extract scalar field data associated with each point on the surface.

        Args:
            fieldname (str): Must be one of: `bipolar_voltage`, `unipolar_voltage`,
                `local_activation_time`, `impedance`, `force`.
            copy (bool, optional):  If True, a copy of the data will be returned. The
                default is False, in which case a view of the data is returned.

        Returns:
            field (np.ndarray): scalar field data
        """

        field = self.fields[fieldname]

        if copy:
            field = field.copy()

        return field

    def add_unipolar_electrograms(
        self,
        unipolar,
        add_bipolar=True,
        add_annotations=True,
        ):
        """Add unipolar electrograms into the Case object.

        The unipolar electrograms and associated data (voltages, names, points)
        are all modified in-place.

        Bipolar electrograms and associated data can also be calculated and modified in-place.

        Args:
            unipolar (np.ndarray): Unipolar electrograms of shape (NxM), where N is the number
                of electrograms and M is the number of time points in each electrogram.
            add_bipolar (bool): If True, bipolar electrograms and associated data
                (voltages, names, points) will be determined from the unipolar
                electrograms and returned.
            add_annotations (bool): If True, default annotations of the electrograms will
                be created and returned. The window of interest will be set to
                cover the entire period of the electrogram traces, and the reference
                annotations will all be set to 0 ms.
        """

        if len(unipolar) != len(self.points):
            raise ValueError(
                "There must be one electrogram per point in the surface. "
                "However, the number of electrograms in the data file is ", len(unipolar),
                " but the number of points in the mesh is ", len(self.points), " . "
            )

        names = np.arange(len(unipolar)).astype(str)
        bipolar, pair_indices = bipolar_from_unipolar_surface_points(
            unipolar=unipolar,
            indices=self.indices,
        )

        pairs_A, pairs_B = pair_indices.T
        unipolar_A = unipolar[pairs_A]
        unipolar_B = unipolar[pairs_B]
        points_A = self.points[pairs_A]
        points_B = self.points[pairs_B]

        # Use both unipolar_A and unipolar_B for `egm`` to mirror the data structure obtained from the clinical cases.
        unipolar_egm = Electrogram(
            egm=np.concatenate([unipolar_A[:, :, np.newaxis], unipolar_B[:, :, np.newaxis]], axis=2),
            points=np.concatenate([points_A[:, :, np.newaxis], points_B[:, :, np.newaxis]], axis=2),
            voltage=np.ptp(unipolar[pair_indices], axis=2),  # axis=2 because of the shape of unipolar[pairs_indices]
            names=np.asarray(['_'.join(pair) for pair in names[pair_indices]])
        )

        self.electric.unipolar_egm = unipolar_egm

        if add_bipolar:

            bipolar_egm = Electrogram(
                egm=bipolar,
                points=self.points,
                voltage=np.ptp(bipolar, axis=1),
                names=names,
            )
            self.electric.bipolar_egm = bipolar_egm

        if add_annotations:
                        
            woi = np.zeros((len(names), 2), dtype=int)
            woi[:, 1] = unipolar.shape[1]

            annotations = Annotations(
                window_of_interest=woi,
                local_activation_time=None,  # TODO: calculate local activation times
                reference_activation_time=np.zeros_like(woi[:, 0], dtype=int)
            )
            self.electric.annotations = annotations
