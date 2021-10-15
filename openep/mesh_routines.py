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
import pyvista
import pymeshfix
import networkx as nx

from matplotlib.cm import jet_r

__all__ = [
    "compute_field",
    "calculate_per_triangle_field",
    "calculate_point_distance_max",
    "calculate_vertex_path",
    "calculate_vertex_distance",
    "calculate_mesh_volume",
    "calculate_field_area",
]


# TODO: remove this function as the `scalars` keyword of pyvista.Plotter().add_mesh can be used instead
def compute_field(
    mesh,
    fieldname,
    minval=0,
    maxval=1,
    color_map=jet_r,
    below_color=(180, 180, 180, 255),
    above_color=(255, 0, 255, 255),
    nan_color=(180, 180, 180, 255),
) -> np.ndarray:
    case = mesh._kwargs["parent_obj"]
    field = case.fields[fieldname]
    new_field = np.zeros((field.shape[0], 4), np.int32)

    for idx in range(len(field)):
        val = field[idx]

        if np.isnan(val):
            col = nan_color
        elif val < minval:
            col = below_color
        elif val > maxval:
            col = above_color
        else:
            valscaled = (val - minval) / (maxval - minval)
            col = color_map(valscaled)

            if isinstance(col[0], float):
                col = [int(c * 255) for c in col]

        new_field[idx] = col

    mesh.visual.vertex_colors[:] = new_field

    return new_field


def repair_mesh(mesh: pyvista.PolyData) -> pyvista.PolyData:
    """
    Repair the given mesh using pymeshfix.
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


def calculate_field_area(
    mesh: pyvista.PolyData, field: np.ndarray, threshold: float
) -> float:
    """
    Calculate the area of triangles whose values are at or below the given threshold.

    Args:
        mesh (PolyData): pyvista mesh
        field ndarray: field used to select triangles
        threshold (float): value at or below which triangles are selected to include in calculation

    Returns:
        float: total area of selected triangles
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

    graph = create_edge_graph(mesh)
    path = nx.shortest_path(graph, source=start_idx, target=end_idx, weight="length")

    return np.array(path, int)


# TODO: This function can be replaced by pyvista.geodesic
#       Just need to handle the exception raised by pyvista.geodesic if there is no path
def create_edge_graph(mesh: pyvista.PolyData) -> nx.Graph:
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
