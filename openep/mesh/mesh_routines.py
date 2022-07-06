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
Analyse a mesh - :mod:`openep.mesh.mesh_routines`
=================================================

This module provides methods for calculating the mesh :ref:`surface area and volume
<geometry>`, calculating :ref:`Euclidian and geodesic distances <distances>`
bewteen points on a mesh, and :ref:`identifying the free boundaries
of a mesh <boundaries>`,

.. _geometry:

Calculating the mesh surface area and volume
--------------------------------------------

.. autofunction:: calculate_mesh_volume

.. autofunction:: calculate_field_area

.. autofunction:: calculate_per_triangle_field

.. _distances:

Distances between points on a mesh
----------------------------------

.. autofunction:: calculate_vertex_distance

.. autofunction:: calculate_vertex_path


.. _boundaries:

Identifying and analysing the free boundaries of a mesh
-------------------------------------------------------

.. autofunction:: get_free_boundaries

.. autoclass:: FreeBoundary
    :members: separate_boundaries, calculate_lengths, calculate_areas

"""

from attr import attrs
from typing import Union

import numpy as np
import scipy.stats

import pyvista
import pymeshfix
import trimesh

__all__ = [
    "get_free_boundaries",
    "calculate_mesh_volume",
    "repair_mesh",
    "calculate_per_triangle_field",
    "calculate_field_area",
    "calculate_vertex_distance",
    "calculate_vertex_path",
    "voxelise",
]


def _create_trimesh(pyvista_mesh):
    """Convert a pyvista mesh into a trimesh mesh.

    Args:
        pyvista_mesh (pyvista.PolyData): The pyvista mesh from which the trimesh mesh will be generated

    Returns:
        trimesh_mesh (trimesh.Trimesh): The generated trimesh mesh
    """

    vertices = pyvista_mesh.points
    faces = pyvista_mesh.faces.reshape(pyvista_mesh.n_faces, 4)[:, 1:]  # ignore to number of vertices per face

    return trimesh.Trimesh(vertices, faces, process=False)


@attrs(auto_attribs=True, auto_detect=True)
class FreeBoundary:
    """
    Class for storing information on the free boundaries of a mesh.

    Args:
        points (np.ndarray): (N,3) array of coordinates
        lines (np.ndarray): (M, 2) array containing the indices of points connected
            by an edge. Values should be in the interval [0, N-1].
        n_boundaries (int): the number of free boundaries
        n_points_per_boundary (np.ndarray): the number of points in each boundary
        original_lines (np.ndarray): (M, 2) array. Same as `lines`, but the indices
            should correspond to those in the original mesh from which free boundaries
            were identified.
    """

    points: np.ndarray
    lines: np.ndarray
    n_boundaries: int
    n_points_per_boundary: np.ndarray
    original_lines: np.ndarray

    def __attrs_post_init__(self):

        if self.n_boundaries == 0:
            self._start_indices = None
            self._stop_indices = None
            self._boundary_meshes = None
            return None

        # We'll use the start and stop indices for separating the (N,2) ndarray
        # of indices into separate arrays for each boundary
        start_indices = list(np.cumsum(self.n_points_per_boundary[:-1] - 1))
        start_indices.insert(0, 0)
        self._start_indices = np.asarray(start_indices)

        stop_indices = start_indices[1:]
        stop_indices.append(None)
        self._stop_indices = np.asarray(stop_indices)

        self._boundary_meshes = None

    def separate_boundaries(self, original_lines=False):
        """
        Creates a list of numpy arrays where each array contains the indices of
        node pairs in a single free boundary.

        Args:
            original_lines (bool):
                If True, `FreeBoundary.original_lines` will be used.
                If False, `FreeBoundary.lines` will be used.

        Returns:
            boundaries (list): a list of numpy arrays - one array per free boundary.
            Each array is of shape Nx2, where N is the number of lines in a given boundary.
            Each array contains the indices of pairs of nodes that make up each line in the
            boundary.

        """

        if self.n_boundaries == 0:
            return np.array([])

        lines = self.original_lines if original_lines else self.lines

        separate_boundaries = [
            lines[start:stop] for start, stop in zip(self._start_indices, self._stop_indices)
        ]

        return separate_boundaries

    def calculate_lengths(self):
        """
        Calculates the length of the perimeter of each free boundary.

        Returns:
            lengths (np.ndarray): the perimeter of each free boundary

        """

        if self.n_boundaries == 0:
            return np.array([])

        lengths = [
            self._line_length(self.points[self.lines[start:stop]]) for
            start, stop in zip(self._start_indices, self._stop_indices)
        ]

        return np.asarray(lengths)

    def _line_length(self, points):
        """
        Calculates the length of a line defined by the positions of a sequence of points.

        Args:
            points (ndarray): Nx3 array of cartesian coordinates of points along the line.

        Returns:
            length (float): length of the line
        """

        distance_between_neighbours = np.sqrt(np.sum(np.square(points[:, 0, :] - points[:, 1, :]), axis=1))
        total_distance = np.sum(distance_between_neighbours)

        return float(total_distance)

    def calculate_areas(self):
        """
        Calculates the cross-sectional area of each boundary.

        Returns:
            areas (np.ndarray): the area of each free boundary

        """

        if self.n_boundaries == 0:
            return np.array([])

        if self._boundary_meshes is None:
            self._create_boundary_meshes()

        areas = [mesh.area for mesh in self._boundary_meshes]

        return np.asarray(areas)

    def _create_boundary_meshes(self):
        """
        Create a pyvista.PolyData mesh for each boundary.

        This determines the geometric centre of the boundary, adds a point at the centre, then
        creates a new mesh in which every edge forms a triangle with the central point. This
        thus creates a surface of the boundary, from which the cross-sectional area can be calculated.
        """

        boundary_meshes = []
        boundaries = self.separate_boundaries(original_lines=False)

        for boundary in boundaries:

            points = self.points[boundary[:, 0]]
            center = np.mean(points, axis=0)
            points = np.vstack([center, points])

            num_points = points.shape[0]
            n_vertices_per_node = np.full(num_points - 1, fill_value=3, dtype=int)

            vertex_one = np.zeros(num_points - 1, dtype=int)  # all triangles include the central point, index 0
            vertex_two = np.arange(1, num_points)
            vertex_three = np.roll(vertex_two, shift=-1)
            faces = np.vstack([n_vertices_per_node, vertex_one, vertex_two, vertex_three]).T.ravel()

            boundary_meshes.append(pyvista.PolyData(points, faces))

        self._boundary_meshes = boundary_meshes

        return None


def get_free_boundaries(mesh):
    """
    Determines the freeboundary/outlines of the 3-D mesh.

    Args:
        mesh (pyvista.PolyData): An open mesh for which the free boundaries will be determined.

    Returns:
        free_boundaries (openep.mesh.FreeBoundary):
            The free boundaries of the open mesh.
    """

    tm_mesh = _create_trimesh(mesh)
    boundaries = tm_mesh.outline().entities

    if boundaries.size == 0:
        return FreeBoundary(
            points=np.array([]),
            lines=np.array([]),
            n_boundaries=0,
            n_points_per_boundary=np.array([]),
            original_lines=np.array([]),
        )

    # determine information about each boundary
    original_indices = np.concatenate([line.points for line in boundaries])
    new_indices = np.arange(original_indices.size)

    n_points_per_boundary = np.asarray([line.points.size for line in boundaries])
    n_boundaries = boundaries.size

    # Create an array pairs of neighbouring nodes for each boundary
    original_lines = np.vstack([original_indices[:-1], original_indices[1:]]).T
    new_lines = np.vstack([new_indices[:-1], new_indices[1:]]).T

    # Ignore the neighbours that are part of different boundaries
    keep_lines = np.full_like(new_lines[:, 0], fill_value=True, dtype=bool)
    keep_lines[n_points_per_boundary[:-1].cumsum()-1] = False

    original_lines = original_lines[keep_lines]
    new_lines = new_lines[keep_lines]

    # Get the {x,y,z} coordinates of the first node in each pair
    points = tm_mesh.vertices[original_indices]

    return FreeBoundary(
        points=points,
        lines=new_lines,
        n_boundaries=n_boundaries,
        n_points_per_boundary=n_points_per_boundary,
        original_lines=original_lines,
    )


def calculate_mesh_volume(
    mesh: pyvista.PolyData,
    fill_holes: bool = True,
) -> float:
    """
    Calculate the volume of a mesh.

    Args:
        mesh (PolyData): mesh for which the volume will be calculated
        fill_holes: if True, holes in the mesh are filled. If holes are present the volume is meaningless unless
        they are filled.

    Returns:
        The volume of the mesh.
    """

    if fill_holes:
        mesh = repair_mesh(mesh)

    return mesh.volume


def repair_mesh(mesh: pyvista.PolyData) -> pyvista.PolyData:
    """
    Fill the holes of a mesh to make it watertight.

    Args:
        mesh (PolyData): mesh to be repaired.

    Returns:
        mesh (PolyData): the repaired mesh.
    """
    mf = pymeshfix.MeshFix(mesh)
    mf.repair()

    return mf.mesh


def calculate_per_triangle_field(mesh: pyvista.PolyData, field: np.ndarray) -> np.ndarray:
    """
    Calculate a per-triangle field from the given per-vertex field. For each triangle the mean of the vertex values is
    calculated as the triangle value.

    Args:
        mesh (PolyData): PolyData mesh
        field: per-vertex field to convert

    Returns:
        np.ndarray per-triangle field
    """
    faces = mesh.faces.reshape(-1, 4)[:, 1:]
    return field[faces].mean(axis=1)


def calculate_field_area(
    mesh: pyvista.PolyData, field: np.ndarray, threshold: float
) -> float:
    """
    Calculate the total surface area of cells whose corresponding values in `field` are
    less than or equal to the given threshold.

    Args:
        mesh (PolyData): pyvista mesh
        field (ndarray): scalar values that will be filtered based on the given threshold
        threshold (float): cells with values in `field` less than or equal to this value
            will be included when calculating the surface area.

    Returns:
        float: total area of selected cells

    Note
    ----
    This function makes use of :func:`openep.mesh.mesh_routines.calculate_per_triangle_field`

    """

    areas = mesh.compute_cell_sizes(
        length=False,
        area=True,
        volume=False,
    )['Area']

    tri_field = calculate_per_triangle_field(mesh, field)
    selection = tri_field <= threshold
    selected_areas = areas[selection]

    return selected_areas.sum()


def calculate_vertex_distance(
    mesh: pyvista.PolyData,
    start_index: int,
    end_index: int,
    metric: str = "geodesic",
) -> float:
    """
    Calculate the distance from vertex at `start_idx` to `end_idx`.

    Either the Euclidian or geodesic distance can be calculated.

    Args:
        mesh (PolyData): Polydata mesh
        start_index (int) : index of starting vertex
        end_index (int) : index of ending vertex
        metric (str): The distance metric to use. The distance function can
            be 'geodesic' or 'euclidian'.

    Returns:
        float: distance between vertices
    """

    if metric not in {"geodesic", "euclidian"}:
        raise ValueError("metric must be on of: geodesic, euclidian")

    if metric == "euclidian":

        distance = np.linalg.norm(
            mesh.points[start_index] - mesh.points[end_index]
        )

        return distance

    try:
        distance = mesh.geodesic_distance(start_index, end_index)
    except(ValueError):
        distance = np.NaN

    return distance


def calculate_vertex_path(
    mesh: pyvista.PolyData,
    start_index: int,
    end_index: int
) -> np.ndarray:
    """
    Calculate the path from vertex at `start_idx` to `end_idx` as a path of vertices through the mesh.

    This is a wrapper around pyvista.PolyData.geodesic, but it returns an empty array if no path
    exist between the two vertices.

    Args:
        mesh (PolyData): Polydata mesh
        start_index (int) : index of starting vertex
        end_index (int) : index of ending vertex

    Returns:
        ndarray: Array of vertex indices defining the path
    """

    try:
        path_mesh = mesh.geodesic(start_vertex=start_index, end_vertex=end_index)
        path = np.asarray(path_mesh.point_data['vtkOriginalPointIds'][path_mesh.lines[1:]])

    except(ValueError):
        path = np.array([])

    return path


def _get_unreferenced_points(mesh):
    """Determine indices of points not referenced in the triangulation"""

    indices = np.arange(mesh.n_points)
    referenced_indices = np.unique(mesh.faces.reshape(mesh.n_faces, 4)[:, 1:].ravel())
    unreferenced_indices = np.isin(indices, referenced_indices, assume_unique=True, invert=True)

    return unreferenced_indices


def _determine_voxel_bins(mesh, edge_length, border=10):
    """Determine the bins to voxelise a mesh.
    
    Args:
        mesh (pyvista.PolyData): Mesh to be voxelised.
        edge_length (float): Edge length of each voxel, in mm.
        border (float): Minimum border around the mesh, in mm. 

    Returns:
        (np.ndarray): Bins that can be used to construct arrays of voxels.

    """

    def _round_up(value, nearest):
        return np.ceil(value / nearest) * nearest

    def _round_down(value, nearest):
        return np.floor(value / nearest) * nearest

    low_values = np.asarray(mesh.bounds[::2])
    high_values = np.asarray(mesh.bounds[1::2])

    low_values = _round_down(low_values, nearest=border) - border / 2
    high_values = _round_up(high_values, nearest=border) + border / 2

    x_low, y_low, z_low = low_values
    x_high, y_high, z_high = high_values

    x_bins = np.arange(x_low, x_high + edge_length, edge_length)
    y_bins = np.arange(y_low, y_high + edge_length, edge_length)
    z_bins = np.arange(z_low, z_high + edge_length, edge_length)

    bin_edges = [x_bins, y_bins, z_bins]
    bin_centres = [
        x_bins[:-1] + edge_length / 2,
        y_bins[:-1] + edge_length / 2,
        z_bins[:-1] + edge_length / 2,
    ]

    return bin_edges, bin_centres


def voxelise(
    mesh: pyvista.PolyData,
    thickness: Union[float, np.ndarray] = 2,
    n_surfaces: int = 11,
    edge_length: float = 1,
    extract_myocardium: bool = False,
) -> pyvista.PolyData:
    """Voxelise a surface mesh.

    Args:
        mesh (PolyData): Surface mesh to be voxelised.
        thickness (float or np.ndarray): If a float, this defines to thickness of the myocardium.
            An array of thicknesses - one per point in the mesh - can be passed to create a voxelised
            mesh with heterogenous thickness.
        edge_length (float): Length of the voxel edges, in mm.
        n_surfaces (int): A series of surface meshes are created by interpolating points between the
            endocardium and epicardium. Points from this series of meshes are used to determine which voxels
            should be filled. The smaller the voxel edge length, the larger :attr:`n_surfaces`
            should be.
        extract_myocardium (bool, optional): If True the voxelised myocardium will be extracted and
            returned. If False, the voxels in a StructuredGrid will be labelled as filled (1) or empty (0),
            and this data stored as point data in the returned mesh.

    Returns:
        StructuredGrid: The voxelised mesh.
    """

    # Don't make any changes to the mesh
    mesh = mesh.copy(deep=True)

    # Remove points not referenced by the triangulation
    not_referenced = _get_unreferenced_points(mesh)
    mesh.remove_points(not_referenced, inplace=True)

    # Compute normals and set thicknesses
    mesh.compute_normals(inplace=True, auto_orient_normals=True, cell_normals=False, point_normals=True)
    thickness = np.full(mesh.n_points, fill_value=thickness, dtype=float) if isinstance(thickness, float) else thickness
    mesh.point_data['Thickness'] = thickness

    # Calculate voxel bins and create output mesh
    bin_edges, bin_centres = _determine_voxel_bins(mesh, edge_length=edge_length)
    bin_centres_x, bin_centres_y, bin_centres_z = bin_centres

    XX, YY, ZZ = np.meshgrid(bin_centres_x, bin_centres_y, bin_centres_z, indexing='xy')  # must use bin centres, must use Cartesian indexing
    voxels = pyvista.StructuredGrid(XX, YY, ZZ)

    voxel_filled = np.zeros(voxels.n_points, dtype=int)  # keep track of which voxels are filled. 0: empty, 1: filled

    # Determine how much to 
    n_voxels_x = np.ptp(bin_centres_x)
    n_voxels_y = np.ptp(bin_centres_y)
    n_voxels_z = np.ptp(bin_centres_z)

    # Create a series - of open meshes and voxelise each mesh
    for shell_distance in np.linspace(0, 1, n_surfaces):

        shell = mesh.copy(deep=True)
        shell.points += mesh.point_data['Normals'] * mesh.point_data['Thickness'][:, np.newaxis] * shell_distance
        shell.subdivide_adaptive(max_edge_len=edge_length, inplace=True)

        mesh_binned = scipy.stats.binned_statistic_dd(
            sample=np.asarray(shell.points),
            values=np.zeros(shell.n_points),
            statistic='count',
            bins=bin_edges,
            expand_binnumbers=True,
        )

        x_indices, y_indices, z_indices = mesh_binned.binnumber - 1
        bin_indices = np.ravel_multi_index(
            [y_indices, x_indices, z_indices],  # we're using Cartesian indexing (as does PyVista)
            dims=np.asarray([n_voxels_y+1, n_voxels_x+1, n_voxels_z+1], dtype=int),
            order='F',
        )

        voxel_filled[bin_indices] = 1

    voxels.point_data['Filled'] = voxel_filled

    if extract_myocardium:
        voxels = voxels.extract_points(voxel_filled.astype(bool))

    return voxels
