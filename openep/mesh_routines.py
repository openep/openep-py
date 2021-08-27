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

from typing import Optional, Union
import numpy as np

from trimesh import Trimesh
import trimesh.repair
import pymeshfix
import networkx as nx

from .case import Case

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


def get_mesh(mesh_case: Union[Case, Trimesh]) -> Trimesh:
    """
    If mesh_case is a Case object, create a mesh from it without backfaces, normals, or recentering. If `mesh_case`
    is a mesh, return this directly. This routine is useful for defining others to accept either type.
    """
    if isinstance(mesh_case, Case):
        return mesh_case.create_mesh(False, False, False)
    else:
        return mesh_case


def repair_mesh(mesh: Trimesh) -> Trimesh:
    """
    Repair the given mesh using pymeshfix.
    """
    mf = pymeshfix.MeshFix(mesh.vertices, mesh.faces)
    mf.repair()

    return Trimesh(mf.v, mf.f)


def create_edge_graph(mesh: Trimesh) -> nx.Graph:
    """
    Create a Graph object of the mesh's edges with the length stored as property 'length'.
    """
    edges = mesh.edges_unique
    edge_lengths = mesh.edges_unique_length

    graph = nx.Graph()
    for edge, length in zip(edges, edge_lengths):
        graph.add_edge(*edge, length=length)

    return graph


def calculate_per_triangle_field(mesh: Trimesh, field: np.ndarray) -> np.ndarray:
    """
    Calculate a per-triangle field from the given per-vertex field. For each triangle the mean of the vertex values is
    calculated as the triangle value.

    Args:
        mesh (obj): Trimesh object
        field: per-vertex field to convert

    Returns:
        np.ndarray per-triangle field
    """
    return field[mesh.faces].mean(axis=1)


def calculate_mesh_volume(
    mesh_case: Union[Case, Trimesh], fill_holes: bool = True
) -> float:
    """
    Calculate the volume of a mesh.

    Args:
        mesh_case (obj): if a Case object, a mesh is created from it, otherwise it is the mesh to calculate the volume for
        fill_holes: if True, holes in the mesh are filled, if holes are present the volume is meaningless

    Returns:
        The volume float value
    """
    mesh = get_mesh(mesh_case)

    # try to fill simple holes, resorting to repair_mesh if needed
    if fill_holes:
        if not trimesh.repair.fill_holes(mesh):
            mesh = repair_mesh(mesh)

    return mesh.volume


def calculate_field_area(
    mesh_case: Union[Case, Trimesh], field: np.ndarray, threshold: float
) -> float:
    """
    Calculate the area of triangles whose values are at or below the given threshold.

    Args:
        mesh_case(obj): Case or Trimesh object
        field: field used to select triangles
        threshold(float): value at or below which triangles are selected to include in calculation

    Returns:
        float: total area of selected triangles
    """
    mesh = get_mesh(mesh_case)
    areas = mesh.area_faces
    tri_field = calculate_per_triangle_field(mesh, field)
    selection = tri_field <= threshold
    selected_areas = areas[selection]

    return selected_areas.sum()


def calculate_vertex_distance(
    mesh_case: Union[Case, Trimesh], start_idx: int, end_idx: int
) -> float:
    """
    Calculate the euclidean distance from vertex at `start_idx` to `end_idx`.

    Args:
        mesh_case(obj): Case or Trimesh object
        start_idx(int): index of starting vertex
        end_idx(int): index of ending vertex

    Returns:
        float: distance between vertices
    """
    mesh = get_mesh(mesh_case)
    start_vertex = mesh.vertices[start_idx]
    end_vertex = mesh.vertices[end_idx]

    return np.linalg.norm(start_vertex - end_vertex)


def calculate_vertex_path(
    mesh_case: Union[Case, Trimesh], start_idx: int, end_idx: int
) -> np.ndarray:
    """
    Calculate the path from vertex at `start_idx` to `end_idx` as a path of vertices through the mesh.

    Args:
        mesh_case(obj): Case or Trimesh object
        start_idx(int): index of starting vertex
        end_idx(int): index of ending vertex

    Returns:
        int: Array of vertex indices defining the path
    """
    mesh = get_mesh(mesh_case)
    graph = create_edge_graph(mesh)

    path = nx.shortest_path(graph, source=start_idx, target=end_idx, weight="length")

    return np.array(path, int)


def calculate_point_distance_max(points, test_points, max_distance):
    results = []

    dists = []

    for p in test_points:
        dist = np.linalg.norm(points - p, axis=1)
        inds = np.argwhere(dist <= max_distance).flatten()
        results.append(inds)
        dists.append(dist)

    return results, dists
