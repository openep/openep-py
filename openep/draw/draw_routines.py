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
import pyvista
import matplotlib.cm

from ..mesh.mesh_routines import FreeBoundary, get_free_boundaries

__all__ = [
    'draw_free_boundaries',
    'draw_map',
]


# TODO: draw_free_boundaries should be an optional parameter to draw_map
#       Make this function private, and call from draw_map
def draw_free_boundaries(
    free_boundaries: FreeBoundary,
    colour: Union[str, List] = "black",
    width: int = 5,
    plotter: pyvista.Plotter = None,
):
    """
    Draw the freeboundaries of a mesh.

    Args:
        free_boundaries (FreeBoundary): FreeBoundary object. Can be generated using
            openep.mesh.free_boundaries.
        colour (str, list): colour or list of colours to render the free boundaries.
        width (int): width of the free boundary lines.
        plotter (pyvista.Plotter): The free boundaries will be added to this plotting object.
            If None, a new plotting object will be created.

    Returns:
        plotter (pyvista.Plotter): Plotting object with the free boundaries added.

    """

    plotter = pyvista.Plotter() if plotter is None else plotter
    colours = [colour] * free_boundaries.n_boundaries if isinstance(colour, str) else colour

    for boundary_index, boundary in enumerate(free_boundaries.separate_boundaries()):

        points = free_boundaries.points[boundary[:, 0]]
        points = np.vstack([points, points[:1]])  # we need to close the loop
        plotter.add_lines(points, color=colours[boundary_index], width=width)

    return plotter


def draw_map(
    mesh: pyvista.PolyData,
    field: np.ndarray,
    plotter: pyvista.Plotter = None,
    add_mesh_kws: dict = None,
    free_boundaries: bool = True,
):
    """
    Project scalar values onto a mesh and optionally draw the free boundaries.

    Args:
        mesh (PolyData): mesh to be drawn
        field (nx1 array): scalar values used to colour the mesh
        plotter (pyvista.Plotter): The mesh will be added to this plotting object.
            If None, a new plotting object will be created.
        add_mesh_kws (dict): Keyword arguments for pyvista.Plotter.add_mesh()
        free_boundaries (bool): If True, the free boundaries will be added to the plot.

    Returns:
        plotter (pyvista.Plotter): Plotting object with the mesh added.
    """

    plotter = pyvista.Plotter() if plotter is None else plotter

    # Create default settings for the plot
    scalar_bar_args = dict(
        interactive=True,
        n_labels=2,
        label_font_size=30,
        below_label=" ",
        above_label=" ",
    )
    if add_mesh_kws is not None and "scalar_bar_args" in add_mesh_kws:
        add_mesh_kws["scalar_bar_args"] = {**scalar_bar_args, **add_mesh_kws["scalar_bar_args"]}

    default_add_mesh_kws = {
        "style": "surface",
        "show_edges": False,
        "smooth_shading": True,
        "annotations": False,
        "cmap": matplotlib.cm.jet_r,
        "clim": (0, 2),
        "above_color": "magenta",
        "below_color": "brown",
        "nan_color": "gray",
        "scalar_bar_args": scalar_bar_args,
    }

    # combine the default and user-given kwargs
    add_mesh_kws = default_add_mesh_kws if add_mesh_kws is None else {**default_add_mesh_kws, **add_mesh_kws}

    plotter.add_mesh(
        mesh=mesh,
        scalars=field,
        **add_mesh_kws,
    )

    if free_boundaries:
        draw_free_boundaries(
            get_free_boundaries(mesh),
            plotter=plotter
        )

    return plotter
