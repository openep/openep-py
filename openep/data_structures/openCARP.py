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

from .._exceptions import NoDataError
from .electric import Electric, Electrogram, Annotations

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

    def __attrs_post_init__(self):

        self.electric = None

    def __repr__(self):
        return f"openCARP mesh with {len(self.points)} points."

    def add_data(self, data_type: str, data_file: str):
        """Add data files to the CARPData.

        Args:
            data_type (str): data type to be added. Must be one of:
                'unipolar'.
            data_file (str): Name of the .dat file containing the data to load.
        """

        supported_types = {
            'unipolar egm',
        }
        if data_type not in supported_types:
            raise ValueError("Unsupported data type: ", data_type)

        func = getattr(self, f"_add_{'_'.join(data_type.split())}_data")
        func(data_file)

    def _add_unipolar_egm_data(self, data_file):

        unipolar = np.loadtxt(data_file)

        if len(unipolar) != len(self.points):
            raise ValueError(
                "The number of points in the data file is ", len(unipolar),
                " but the number of points in the mesh is ", len(self.points), " ."
            )

        if self.electric is not None:
            raise NotImplementedError("Adding multiple data files is not yet supported.")

        # determine the bipolar neibhours
        names = np.arange(len(unipolar)).astype(str)
        bipolar, pair_indices = self.bipolar_from_unipolar(unipolar)
        
        bipolar_egm = Electrogram(
            egm=bipolar,
            points=self.points,
            voltage=np.ptp(bipolar, axis=1),
            names=names,
        )

        unipolar_egm = Electrogram(
            egm=unipolar[pair_indices],
            points=self.points[pair_indices],
            voltage=np.ptp(unipolar[pair_indices], axis=1),
            names=np.asarray(['_'.join(pair) for pair in names[pair_indices]])
        )

        woi = np.zeros((len(names), 2), dtype=int)
        woi[:, 1] = unipolar.shape[1]
        annotations = Annotations(
            window_of_interest=woi,
            local_activation_time=None,  # TODO: calculate local activation times
            reference_activation_time=np.zeros_like(woi[:, 0], dtype=int)
        )

        self.electric = Electric(
            names=names,
            internal_names=names,
            bipolar_egm=bipolar_egm,
            unipolar_egm=unipolar_egm,
            reference_egm=None,
            ecg=None,
            impedance=None,
            surface=None,
            annotations=annotations,
        )

        self._unipolar = unipolar
        self._unipolar_pairs = pair_indices

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

    def bipolar_from_unipolar(self, unipolar):
        """
        Calculate bipolar electrograms from unipolar electrograms for each point on the mesh.

        Args:
            unipolar (np.ndarray): Unipolar electrograms

        Returns
            bipolar (np.ndarray): Bipolar electrograms
        """

        bipolar = np.full_like(unipolar, fill_value=np.NaN)
        pair_indices = np.full((len(unipolar), 2), fill_value=0, dtype=int)

        for index, index_unipolar in enumerate(unipolar):

            connected_vertices = self._find_connected_vertices(self.indices, index)
            index_bipolar, pair_index = self._bipolar_from_unipolar(
                unipolar=index_unipolar,
                neighbours=unipolar[connected_vertices],
            )
            bipolar[index] = index_bipolar
            pair_indices[index] = [index, pair_index]

        return bipolar, pair_indices

    def _find_connected_vertices(self, faces, index):
        """
        Find all points connected to a given point by a single edge.

        Adapted from: https://github.com/pyvista/pyvista-support/issues/96#issuecomment-571864471

        Args:
            faces (np.ndarray): faces of a mesh
            index (int): index of point for which we want to find the neighbouring points
        """

        connected_faces = [i for i, face in enumerate(faces) if index in face]
        connected_vertices = np.unique(faces[connected_faces])

        return connected_vertices[connected_vertices != index]

    def _bipolar_from_unipolar(self, unipolar, neighbours):
        """
        Calculate the bipolar electrogram of a given point from a series of unipolar electrograms.

        Args:
            unipolar (np.ndarray): unipolar electrogram at a given point
            neighbours (np.narray): unipolar electrograms at all points
                neighbouring the given point
        """

        difference = neighbours - unipolar
        voltage = np.ptp(difference, axis=1)
        pair_index = np.argmax(voltage)
        bipolar_electrogram = difference[pair_index]

        return bipolar_electrogram, pair_index
