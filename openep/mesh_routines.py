from typing import Optional, Union
import numpy as np

from trimesh import Trimesh
import trimesh.repair
import pymeshfix
import networkx as nx

from .case import Case

from matplotlib.cm import jet_r


def compute_field(
        mesh,
        fieldname,
        minval=0,
        maxval=1,
        color_map=jet_r,
        below_color=(0, 0, 0, 255),
        above_color=(255, 0, 255, 255),
        nan_color=(50, 50, 50, 255),
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
    If `mesh_case` is a Case object, create a mesh from it without backfaces, normals, or recentering. If `mesh_case`
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
        mesh: Trimesh object
        field: per-vertex field to convert

    Returns:
        np.ndarray per-triangle field
    """
    return field[mesh.faces].mean(axis=1)


def calculate_mesh_volume(mesh_case: Union[Case, Trimesh], fill_holes: bool = True) -> float:
    """
    Calculate the volume of a mesh.

    Args:
        mesh_case: if a Case object, a mesh is created from it, otherwise it is the mesh to calculate the volume for
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


def calculate_field_area(mesh_case: Union[Case, Trimesh], field: np.ndarray, threshold: float) -> float:
    """
    Calculate the area of triangles whose values are at or below the given threshold.

    Args:
        mesh_case: Case or Trimesh object
        field: field used to select triangles
        threshold: value at or below which triangles are selected to include in calculation

    Returns:
        Float total area of selected triangles
    """
    mesh = get_mesh(mesh_case)
    areas = mesh.area_faces
    tri_field = calculate_per_triangle_field(mesh, field)
    selection = tri_field <= threshold
    selected_areas = areas[selection]

    return selected_areas.sum()


def calculate_vertex_distance(mesh_case: Union[Case, Trimesh], start_idx: int, end_idx: int) -> float:
    """
    Calculate the euclidean distance from vertex at `start_idx` to `end_idx`.

    Args:
        mesh_case: Case or Trimesh object
        start_idx: index of starting vertex
        end_idx: index of ending vertex

    Returns:
        Float distance between vertices
    """
    mesh = get_mesh(mesh_case)
    start_vertex = mesh.vertices[start_idx]
    end_vertex = mesh.vertices[end_idx]

    return np.linalg.norm(start_vertex - end_vertex)


def calculate_vertex_path(mesh_case: Union[Case, Trimesh], start_idx: int, end_idx: int) -> np.ndarray:
    """
    Calculate the path from vertex at `start_idx` to `end_idx` as a path of vertices through the mesh.

    Args:
        mesh_case: Case or Trimesh object
        start_idx: index of starting vertex
        end_idx: index of ending vertex

    Returns:
        Array of vertex indices defining the path
    """
    mesh = get_mesh(mesh_case)
    graph = create_edge_graph(mesh)

    path = nx.shortest_path(graph, source=start_idx, target=end_idx, weight='length')

    return np.array(path, int)
