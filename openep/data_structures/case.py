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


from attr import attrs
from typing import Any, Dict, Optional, Tuple, List

import numpy as np
import pyvista

from .surface import Fields
from .electric import Electric

__all__ = []


@attrs(auto_attribs=True, auto_detect=True)
class Case:
    """OpenEP Case object."""
    name: str
    points: np.ndarray
    indices: np.ndarray
    fields: Fields
    electric: Electric
    rf: Optional[Dict[str, Any]] = None
    notes: Optional[List] = None

    def __repr__(self):
        return f"{self.name}( nodes: {self.points.shape} indices: {self.indices.shape} {self.fields} )"

    def create_mesh(
        self, vertex_norms: bool = False,
        recenter: bool = True,
        back_faces: bool = False,
    ):
        """
        Create a new mesh object from the stored nodes and indices

        Args:
            vertex_norms: if True, calculate vertex normals in the mesh
            recenter: if True, recenter the mesh on the origin
            back_faces: if True, calculate back face triangles
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

        # TODO: Refactor to remove this statement and the vertex_norms parameter.
        if vertex_norms:
            _ = mesh.point_normals  # compute vertex normals

        if recenter:
            mesh.translate(
                -np.asarray(mesh.center)
            )  # recenter mesh to origin, helps with lighting in default scene

        return mesh

    def get_surface_data(self, copy: bool = False) -> Tuple[np.ndarray, np.ndarray]:
        """
        Returns the node and triangle index matrices, copying them if `copy` is True.
        """
        nodes = self.points
        indices = self.indices

        if copy:
            nodes = nodes.copy()
            indices = indices.copy()

        return nodes, indices

    def get_field(self, fieldname: str, copy: bool = False) -> np.ndarray:
        """
        Returns the named field array, copying if `copy` is True.
        """

        field = self.fields[fieldname]

        if copy:
            field = field.copy()

        return field
