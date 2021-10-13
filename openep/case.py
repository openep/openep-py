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


import numpy as np
from typing import Any, Dict, Optional, Tuple
import pyvista

__all__ = ["Case"]


class Case:
    def __init__(
        self,
        name: str,
        nodes: np.ndarray,
        indices: np.ndarray,
        fields: Dict[str, np.ndarray],
        electric: Dict[str, np.ndarray],
        rf: Dict[str, np.ndarray],
        other_data: Optional[Dict[str, Any]] = None,
    ):
        self.name: str = name
        self.nodes: np.ndarray = nodes
        self.indices: np.ndarray = indices
        self.fields: Dict[str, np.ndarray] = dict(fields)
        self.electric: Dict[str, np.ndarray] = dict(electric)
        self.rf: Dict[str, np.ndarray] = dict(rf)
        self.other_data: Dict[str, Any] = dict(other_data or {})

    def __repr__(self):
        return f"{self.name}( nodes: {self.nodes.shape} indices: {self.indices.shape} fields: {tuple(self.fields)} )"

    def create_mesh(
        self, vertex_norms: bool = True, recenter: bool = True, back_faces: bool = True
    ):
        """
        Create a new mesh object from the stored nodes and indices

        Args:
            vertex_norms: if True, calculate vertex normals in the mesh
            recenter: if True, recenter the mesh on the origin
            back_faces: if True, calculate back face triangles
        """
        if back_faces:
            
            inds = self.indices.reshape(-1, 4)            
            inds_inverted = inds[:, [0, 1, 3, 2]]
            
            inds = np.concatenate(
                [inds.ravel(), inds_inverted.ravel()]
            )  # include each triangle twice for both surfaces
        else:
            inds = self.indices

        mesh = pyvista.PolyData(self.nodes, inds)
        
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
        nodes = self.nodes
        indices = self.indices

        if copy:
            nodes = nodes.copy()
            indices = indices.copy()

        return nodes, indices

    def get_field(self, fieldname: str, copy: bool = False) -> np.ndarray:
        """
        Returns the named field array, copying if `copy` is True.
        """
        if fieldname not in self.fields:
            raise ValueError(f"Field '{fieldname}' not found")

        field = self.fields[fieldname]

        if copy:
            field = field.copy()

        return field
