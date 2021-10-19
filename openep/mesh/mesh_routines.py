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

from dataclasses import dataclass

import numpy as np
import networkx as nx

import pyvista
import pymeshfix
import trimesh

__all__ = [
    "calculate_per_triangle_field",
    "calculate_mesh_volume",
    "calculate_field_area",
    "calculate_vertex_distance",
    "calculate_vertex_path",
    "calculate_point_distance_max",
    "get_free_boundaries",
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


@dataclass
class FreeBoundary:
    """Class for storing information on the free boundaries of a mesh."""

    points: np.ndarray
    lines: np.ndarray
    n_boundaries: int
    n_points_per_boundary: np.ndarray
    original_lines: np.ndarray

    def __post_init__(self):

        start_indices = list(np.cumsum(self.n_points_per_boundary[:-1] - 1))
        start_indices.insert(0, 0)
        self._start_indices = np.asarray(start_indices)

        stop_indices = start_indices[1:]
        stop_indices.append(None)
        self._stop_indices = np.asarray(stop_indices)

        self._boundary_meshes = None

    def separate_boundaries(self, original_lines=False):
        """
        Returns a list of numpy arrays where each array contains the indices of
        node pairs in a single free boundary.

        Args:
            original_lines (bool):
                If True, FreeBoundary.original_indices will be used.
                If False, FreeBoundary.original_indices will be used.

        """

        lines = self.original_lines if original_lines else self.lines

        separate_boundaries = [
            lines[start:stop] for start, stop in zip(self._start_indices, self._stop_indices)
        ]

        return separate_boundaries

    def calculate_lengths(self):
        """
        Returns the length of the perimeter of each free boundary.
        """

        lengths = [
            self._line_length(self.points[start:stop]) for
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

        distance_between_neighbours = np.sqrt(np.sum(np.square(np.diff(points, axis=0)), axis=1))
        total_distance = np.sum(distance_between_neighbours)

        return float(total_distance)

    def calculate_areas(self):
        """
        Returns the total area of the faces in each boundary.
        """

        if self._boundary_meshes is None:
            self._create_boundary_meshes()

        areas = [mesh.area for mesh in self._boundary_meshes]

        return np.asarray(areas)

    def _create_boundary_meshes(self):
        """
        Create a pyvista.PolyData mesh for each boundary.
        """

        boundary_meshes = []
        boundaries = self.separate_boundaries(original_lines=False)

        for boundary in boundaries:

            points = self.points[boundary[:, 0]]
            center = np.mean(points, axis=0)
            points = np.vstack([center, points])

            num_points = points.shape[0]
            n_vertices_per_node = np.full(num_points - 1, fill_value=3, dtype=int)
            vertex_one = np.zeros(num_points - 1, dtype=int)  # all triangles include the central point
            vertex_two = np.arange(1, num_points)
            vertex_three = np.roll(vertex_two, shift=-1)
            faces = np.vstack([n_vertices_per_node, vertex_one, vertex_two, vertex_three]).T.ravel()

            boundary_meshes.append(pyvista.PolyData(points, faces))

        self._boundary_meshes = boundary_meshes

        return None


def get_free_boundaries(mesh):
    """Gets the freeboundary/outlines of the 3-D mesh and returns the indices.

    Args:
        mesh (pyvista.PolyData): Open mesh for which the free boundaries will be determined.

    Returns:
        freeboundaries_info (dict):
            * freeboundary, Nx2 numpy array of indices of neighbouring points in the free boundaries
    """

    tm_mesh = _create_trimesh(mesh)

    # extract the boundary information
    boundaries = tm_mesh.outline()
    boundaries_lines = boundaries.entities

    # determine information about each boundary
    original_indices = np.concatenate([line.points for line in boundaries_lines])
    n_points_per_boundary = np.asarray([line.points.size for line in boundaries_lines])
    n_boundaries = tm_mesh.outline().entities.size
    indices = np.arange(original_indices.size)

    # Create an array pairs of neighbouring nodes for each boundary
    original_lines = np.vstack([original_indices[:-1], original_indices[1:]]).T
    lines = np.vstack([indices[:-1], indices[1:]]).T

    # Ignore the neighbours that are part of different boundaries
    keep_lines = np.full_like(lines[:, 0], fill_value=True, dtype=bool)
    keep_lines[n_points_per_boundary[:-1].cumsum()-1] = False
    original_lines = original_lines[keep_lines]
    lines = lines[keep_lines]

    # Get the {x,y,z} coordinates of the first node in each pair
    points = tm_mesh.vertices[original_indices]

    return FreeBoundary(
        points=points,
        lines=lines,
        n_boundaries=n_boundaries,
        n_points_per_boundary=n_points_per_boundary,
        original_lines=original_lines,
    )


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
        mesh = _repair_mesh(mesh)

    return mesh.volume


def _repair_mesh(mesh: pyvista.PolyData) -> pyvista.PolyData:
    """
    Repair the given mesh using pymeshfix.
    """
    mf = pymeshfix.MeshFix(mesh)
    mf.repair()

    return mf.mesh


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
    start_idx: int,
    end_idx: int
) -> float:
    """
    Calculate the euclidean distance from vertex at `start_idx` to `end_idx`.

    Args:
        mesh (PolyData): Polydata mesh
        start_index (int) : index of starting vertex
        end_index (int) : index of ending vertex

    Returns:
        float: distance between vertices
    """
    start_vertex = mesh.points[start_idx]
    end_vertex = mesh.points[end_idx]

    return np.linalg.norm(start_vertex - end_vertex)


def calculate_vertex_path(
    mesh: pyvista.PolyData,
    start_idx: int,
    end_idx: int
) -> np.ndarray:
    """
    Calculate the path from vertex at `start_idx` to `end_idx` as a path of vertices through the mesh.

    Args:
        mesh (PolyData): Polydata mesh
        start_index (int) : index of starting vertex
        end_index (int) : index of ending vertex

    Returns:
        ndarray: Array of vertex indices defining the path
    """

    graph = _create_edge_graph(mesh)
    path = nx.shortest_path(graph, source=start_idx, target=end_idx, weight="length")

    return np.array(path, int)


# TODO: This function can be replaced by pyvista.geodesic
#       Just need to handle the exception raised by pyvista.geodesic if there is no path
def _create_edge_graph(mesh: pyvista.PolyData) -> nx.Graph:
    """
    Create a Graph object of the mesh's edges with the length stored as property 'length'.
    """

    # Avoid creating the graph unnecessarily
    if hasattr(mesh, "_graph"):
        return mesh._graph

    faces = mesh.faces.reshape(-1, 4)[:, 1:]

    edges = np.vstack(
        [
            faces[:, :2],
            faces[:, 1:3],
            faces[:, ::2],
        ]
    )
    edges.sort()  # ensure the lower index is in the first column
    unique_edges = np.unique(edges, axis=0)

    edge_lengths = np.linalg.norm(mesh.points[unique_edges[:, 0]] - mesh.points[unique_edges[:, 1]], axis=1)

    graph = nx.Graph()
    graph.add_nodes_from(np.arange(mesh.n_points))  # ensure all nodes are present
    for edge, length in zip(unique_edges, edge_lengths):
        graph.add_edge(*edge, length=length)

    mesh._graph = graph

    return graph


def calculate_point_distance_max(points, test_points, max_distance):
    results = []

    dists = []

    for p in test_points:
        dist = np.linalg.norm(points - p, axis=1)
        inds = np.argwhere(dist <= max_distance).flatten()
        results.append(inds)
        dists.append(dist)

    return results, dists
