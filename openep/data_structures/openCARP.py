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

"""Module containing classes for storing data from openCARP simulations."""

from attr import attrs
import pyvista
import numpy as np

__all__ = []


@attrs(auto_attribs=True, auto_detect=True)
class CARPData:
    """
    Class for storing information about meshes from openCARP simulations.

    Args:
        points (np.ndarray): 3D coordinates of points on the mesh.
        indices (ndarray): Indices of points that make up each face of the mesh
        unipolar_egm (np.ndarray, optional): Unipolar electrograms for each point.
    """

    points: np.ndarray
    indices: np.ndarray
    unipolar_egm: np.ndarray = None

    def __attrs_post_init__(self):

        self.bipolar_egm = self.bipolar_from_unipolar()
        self.bipolar_voltage = np.ptp(self.bipolar_egm, axis=1)
        self.unipolar_voltage = np.ptp(self.unipolar_egm, axis=1)

    def __repr__(self):
        return f"openCARP mesh with {len(self.unipolar_egm)} points."

    def create_mesh(
        self,
        recenter: bool = True,
    ):
        """
        Create a new mesh object from the stored nodes and indices

        Args:
            recenter: if True, recenter the mesh to the origin

        Returns:
            mesh (pyvista.Polydata): a mesh created from the points and indices
        """

        indices = self.indices
        num_points_per_face = np.full(shape=(len(indices)), fill_value=3, dtype=int)  # all faces have three vertices
        faces = np.concatenate([num_points_per_face[:, np.newaxis], indices], axis=1)

        mesh = pyvista.PolyData(self.points, faces.ravel())

        if recenter:
            mesh.translate(
                -np.asarray(mesh.center)
            )  # recenter mesh to origin, helps with lighting in default scene

        return mesh

    def bipolar_from_unipolar(self):
        """
        Calculate bipolar electrograms from unipolar electrograms for each point on the mesh.
        
        Returns
            all_bipolar_egm (np.ndarray): Bipolar electrograms
        """
        
        all_bipolar_egm = np.full_like(self.unipolar_egm, fill_value=np.NaN)

        for index, egm in enumerate(self.unipolar_egm):

            connected_vertices = self._find_connected_vertices(self.indices, index)
            bipolar_egm = self._bipolar_from_unipolar(
                egm=egm,
                neighbour_egms=self.unipolar_egm[connected_vertices],
            )
            all_bipolar_egm[index] = bipolar_egm

        return all_bipolar_egm

    def _find_connected_vertices(self, faces, index):
        """
        Find all points connected to a given point by a single edge.

        Args:
            faces (np.ndarray): faces of a mesh
            index (int): index of point for which we want to find the neighbouring points
        """

        connected_faces = [i for i, face in enumerate(faces) if index in face]
        connected_vertices = np.unique(faces[connected_faces])

        return connected_vertices[connected_vertices != index]

    def _bipolar_from_unipolar(self, egm, neighbour_egms):
        """
        Calculate the bipolar electrogram of a given point from a series of unipolar electrograms.

        Args:
            egm (np.ndarray): unipolar electrogram at a given point
            neighbour_egms (np.narray): unipolar electrograms at all points
                neighbouring the given point
        """
        
        difference = neighbour_egms - egm
        voltage = np.ptp(difference, axis=1)
        pair_index = np.argmax(voltage)
        bipolar_electrogram = difference[pair_index]

        return bipolar_electrogram
