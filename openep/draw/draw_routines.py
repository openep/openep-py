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
from typing import List, Union

import numpy as np
import trimesh as tm
import pyvista as pv

__all__ = [
    'get_freeboundaries',
    'draw_free_boundaries',
    'draw_map',
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

    return tm.Trimesh(vertices, faces, process=False)


def _create_pymesh(trimesh_mesh):
    """Convert a pyvista mesh into a trimesh mesh.

    Args:
        pyvista_mesh (pyvista.PolyData): The pyvista mesh from which the trimesh mesh will be generated

    Returns:
        trimesh_mesh (trimesh.Trimesh): The generated trimesh mesh
    """

    return pv.wrap(trimesh_mesh)


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

            boundary_meshes.append(pv.PolyData(points, faces))

        self._boundary_meshes = boundary_meshes

        return None


def get_freeboundaries(mesh):
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


# TODO: draw_free_boundaries should be an optional parameter to draw_map
#       Make this function private, and call from draw_map
def draw_free_boundaries(
    free_boundaries: FreeBoundary,
    colour: Union[str, List] = "black",
    width: int = 10,
    plotter: pv.Plotter = None,
):
    """
    Draw the freeboundaries of a mesh.

    Args:
        free_boundaries (FreeBoundary): FreeBoundary object. Can be generated using
            openep.draw_routines.get_free_boundaries.
        colour (str, list): colour or list of colours to render the free boundaries.
        width (int): width of the free boundary lines.
        plotter (pyvista.Plotter): The free boundaries will be added to this plotting object.
            If None, a new plotting object will be created.

    Returns:
        plotter (pyvista.Plotter): Plotting object with the free boundaries added.

    """

    plotter = pv.Plotter() if plotter is None else plotter
    colours = [colour] * free_boundaries.n_boundaries if isinstance(colour, str) else colour

    for boundary_index, boundary in enumerate(free_boundaries.separate_boundaries()):

        points = free_boundaries.points[boundary[:, 0]]
        points = np.vstack([points, points[:1]])  # we need to close the loop
        plotter.add_lines(points, color=colours[boundary_index], width=width)

    return plotter


# TODO: draw_free_boundaries should be a keyword argument
# TODO: should take a pyvista.Plotter object as an optional argument
def draw_map(
    mesh,
    volt,
    freeboundary_color,
    cmap,
    freeboundary_width,
    minval,
    maxval,
    volt_below_color,
    volt_above_color,
    nan_color,
    plot,
    plotter=None,
    **kwargs
):
    """
    plots an OpenEp Voltage Map
    Args:
        mesh (PolyData): mesh to be drawn
        volt (nx1 array): interpolated voltagae values.
        freeboundary_color(str or rgb list): color of the freeboundaries.
        cmap (str): name of the colormap, for eg: jet_r.
        freeboundary_width (float): width of the freeboundary line.
        minval(float): Voltage lower threshold value.
        maxval(float): Voltage upper threshold value.
        volt_below_color(str or 3 item list): Color for all the voltage values below lower threshold.
        volt_above_color(str or 3 item list): Color for all the voltage values above upper threshold.
        nan_color(str or 3 item list): Color for all the nan voltage values in the openep dataset.
        plot (boolean): True to plot, False otherwise.
        **kwargs: Arbitrary keyword arguments.

    Returns:
        obj: p, VTK actor of the mesh.
        obj: pyvista-mesh, Pyvista PolyData object, triangulated surface object from Numpy arrays of the vertices and faces.
        str or nx1 array: volt, 'bip' or interpolated voltagae values.
        str or 3 item list: nan_color, Color for all the nan voltage values in the openep dataset.
        float: minval,Voltage lower threshold value.
        float: maxval,Voltage upper threshold value.
        str: cmap,name of the colormap.
        str or 3 item list: volt_below_color,Color for all the voltage values below lower threshold
        str or 3 item list: volt_above_color,Color for all the voltage values above upper threshold
    """

    volt = mesh.fields[volt] if isinstance(volt, str) else volt
    plotter = pv.Plotter() if plotter is None else plotter
    
    # Plot OpenEp mesh
    sargs = dict(
        interactive=True,
        n_labels=2,
        label_font_size=18,
        below_label="  ",
        above_label="  ",
    )

    freeboundaries = get_freeboundaries(mesh)
    plotter.add_mesh(
        mesh,
        scalar_bar_args=sargs,
        show_edges=False,
        smooth_shading=True,
        scalars=volt,
        nan_color=nan_color,
        clim=[minval, maxval],
        cmap=cmap,
        below_color=volt_below_color,
        above_color=volt_above_color,
    )

    draw_free_boundaries(
        freeboundaries,
        colour=freeboundary_color,
        width=freeboundary_width,
        plotter=plotter
    )

    if plot:
        plotter.show()

    return {
        "hsurf": plotter,
        "pyvista-mesh": mesh,
        "volt": volt,
        "nan_color": nan_color,
        "minval": minval,
        "maxval": maxval,
        "cmap": cmap,
        "volt_below_color": volt_below_color,
        "volt_above_color": volt_above_color,
    }
