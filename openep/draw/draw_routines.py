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

from typing import List, Union

import numpy as np
import pyvista as pv

from ..mesh.mesh_routines import (FreeBoundary, free_boundaries)

__all__ = [
    'draw_free_boundaries',
    'draw_map',
]


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

    freeboundaries = free_boundaries(mesh)
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
