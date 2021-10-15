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

from openep import case_routines

# TOD0: remove global variables
# GLOBAL VARIABLES
plot = False
volt = 0

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
    
    mesh: pv.PolyData
    freeboundary: np.ndarray
    free_boundary_points: np.ndarray
    n_boundaries: int
    n_nodes_per_boundary: np.ndarray
    
    def __post_init__(self):
        
        start_indices = list(np.cumsum(self.n_nodes_per_boundary[:-1] - 1))
        start_indices.insert(0, 0)
        self._start_indices = np.asarray(start_indices)

        stop_indices = start_indices[1:]
        stop_indices.append(None)
        self._stop_indices = np.asarray(stop_indices)
        
        self._boundary_meshes = None
    
    def separate_boundaries(self):
        """
        Returns a list of numpy arrays where each array contains the indices of
        node pairs in a single free boundary.
        
        """
        
        separate_boundaries = [
            self.freeboundary[start:stop] for start, stop in zip(self._start_indices, self._stop_indices)
        ]
        
        return separate_boundaries
        

    def calculate_lengths(self):
        """
        Returns the length of the perimeter of each free boundary.
        """
        
        lengths = [
            self._line_length(self.free_boundary_points[start:stop]) for
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
        boundaries = self.separate_boundaries()
        
        for boundary in boundaries:
            
            points = self.mesh.points[boundary[:, 0]]
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
    all_boundaries_nodes = np.concatenate([line.points for line in boundaries_lines])
    n_nodes_per_boundary = np.asarray([line.points.size for line in boundaries_lines])
    n_boundaries = tm_mesh.outline().entities.size

    # Create an array pairs of neighbouring nodes for each boundary
    free_boundaries = np.vstack([all_boundaries_nodes[:-1], all_boundaries_nodes[1:]]).T

    # Ignore the neighbours that are part of different boundaries
    keep_neighbours = np.full_like(free_boundaries[:, 0], fill_value=True, dtype=bool)
    keep_neighbours[n_nodes_per_boundary[:-1].cumsum()-1] = False
    free_boundaries = free_boundaries[keep_neighbours]
    
    # Get the {x,y,z} coordinates of the first node in each pair
    free_boundaries_points = tm_mesh.vertices[free_boundaries[:, 0]]
    
    # TODO: return a FreeBoundary object rather than a dictionary
    return {
        "freeboundary": free_boundaries,
        "free_boundary_points": free_boundaries_points,
        "n_boundaries": n_boundaries,
        "n_nodes_per_boundary": n_nodes_per_boundary,
    }

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
            If None, a new plotting object will be created and the mesh associated with the
            free boundaries will also be plotted.

    Returns:
        plotter (pyvista.Plotter): Plotting object with the free boundaries added.
        
    """
    
    if plotter is None:
        
        plotter = pv.Plotter()
        plotter.add_mesh(
            free_boundaries.mesh,
            color=colour,
            smooth_shading=True,
            show_edges=False,
            use_transparency=False,
            opacity=0.2,
            lighting=False
        )
    
    colours = [colour] * 7 if isinstance(colour, str) else colour  # there are 7 anatomical regions
    for boundary_index, boundary in enumerate(free_boundaries.separate_boundaries()):
        
        points = free_boundaries.mesh.points[boundary[:, 0]]
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
    **kwargs
):
    """
    plots an OpenEp Voltage Map
    Args:
        mesh (PolyData): mesh to be drawn
        volt (str or nx1 array): 'bip' or interpolated voltagae values.
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

    if volt == "bip":
        volt = mesh.fields["bip"]
    else:
        volt = volt

    # pyvista plotter
    p = pv.Plotter()

    if plot:

        # Plot OpenEp mesh
        sargs = dict(
            interactive=True,
            n_labels=2,
            label_font_size=18,
            below_label="  ",
            above_label="  ",
        )

        freebound = get_freeboundaries(mesh)
        p.add_mesh(
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
        
        freeboundaries =  FreeBoundary(mesh, **freebound)
        draw_free_boundaries(
            freeboundaries,
            colour=freeboundary_color,
            width=freeboundary_width,
            plotter=p
        )
        
        p.show()

    return {
        "hsurf": p,
        "pyvista-mesh": mesh,
        "volt": volt,
        "nan_color": nan_color,
        "minval": minval,
        "maxval": maxval,
        "cmap": cmap,
        "volt_below_color": volt_below_color,
        "volt_above_color": volt_above_color,
    }


def get_voltage_electroanatomic(mesh_case):
    distance_thresh = 10
    rbf_constant_value = 1

    # Anatomic descriptions (Mesh) - nodes and indices
    pts = mesh_case.nodes

    # Electric data
    # Locations – Cartesian co-ordinates, projected on to the surface
    locations = case_routines.get_electrogram_coordinates(mesh_case, "type", "bip")

    i_egm = mesh_case.electric["egm"].T
    i_vp = case_routines.get_mapping_points_within_woi(mesh_case)
    # macthing the shape of ivp with data
    i_vp_egm = np.repeat(i_vp, repeats=i_egm.shape[1], axis=1)
    # macthing the shape of ivp with coords
    i_vp_locations = np.repeat(i_vp, repeats=locations.shape[1], axis=1)

    # Replacing the values outside the window of interest with Nan values
    i_egm[~i_vp_egm] = np.nan
    locations[~i_vp_locations] = np.nan

    # For each mapping point, n, find the voltage amplitude
    max_volt = np.amax(a=i_egm, axis=1).reshape(len(i_egm), 1)
    min_volt = np.amin(a=i_egm, axis=1).reshape(len(i_egm), 1)

    amplitude_volt = np.subtract(max_volt, min_volt)

    for indx in range(amplitude_volt.shape[1]):
        temp_data = amplitude_volt[:, indx]
        temp_coords = locations
        i_nan = np.isnan(temp_data)
        temp_data = temp_data[~i_nan]
        temp_coords = temp_coords[~i_nan]

        interp = case_routines.OpenEPDataInterpolator(
            method="rbf",
            distanceThreshold=distance_thresh,
            rbfConstant=rbf_constant_value,
        )
        vertex_voltage_data = interp.interpolate(x0=temp_coords, d0=temp_data, x1=pts)

    return vertex_voltage_data
