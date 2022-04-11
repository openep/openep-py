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
Draw meshes and plot electrical data - :mod:`openep.draw.draw_routines`
=======================================================================

This module provides functions for :ref:`drawing 3D meshes <mesh>` using `pyvista`
and :ref:`plotting electrical data <electrical>` using `matplotlib.`

.. _mesh:

Draw 3D meshes
--------------

.. autofunction:: draw_map

.. autofunction:: draw_free_boundaries

.. _electrical:

Plotting electrical data
------------------------

.. autofunction:: plot_electrograms

"""

from typing import List, Union

import numpy as np
import pyvista
import matplotlib.cm
import matplotlib.pyplot as plt

from ..mesh.mesh_routines import FreeBoundary, get_free_boundaries

__all__ = [
    'draw_free_boundaries',
    'draw_map',
    'plot_electrograms',
]


def draw_free_boundaries(
    free_boundaries: FreeBoundary,
    colour: Union[str, List] = "black",
    width: int = 5,
    plotter: pyvista.Plotter = None,
    names: List[str] = None,
):
    """
    Draw the freeboundaries of a mesh.

    Args:
        free_boundaries (FreeBoundary): `FreeBoundary` object. Can be generated using
            :func:`openep.mesh.mesh_routines.get_free_boundaries`.
        colour (str, list): colour or list of colours to render the free boundaries.
        width (int): width of the free boundary lines.
        plotter (pyvista.Plotter): The free boundaries will be added to this plotting object.
            If None, a new plotting object will be created.
        names (List(str)): List of names to associated with the actors. Default is None, in which
            case actors will be called 'free_boundary_n', where n is the index of the boundary.

    Returns:
        plotter (pyvista.Plotter): Plotting object with the free boundaries added.

    """

    plotter = pyvista.Plotter() if plotter is None else plotter
    colours = [colour] * free_boundaries.n_boundaries if isinstance(colour, str) else colour
    names = [f"free_boundary_{boundary_index:d}" for boundary_index in range(free_boundaries.n_boundaries)]

    for boundary_index, boundary in enumerate(free_boundaries.separate_boundaries()):

        points = free_boundaries.points[boundary[:, 0]]
        points = np.vstack([points, points[:1]])  # we need to close the loop
        plotter.add_lines(
            points,
            color=colours[boundary_index],
            width=width,
            name=names[boundary_index],
        )

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
    default_scalar_bar_args = dict(
        interactive=False,
        color="#363737",  # set the colour of the text
        title_font_size=12,
        label_font_size=11,
        n_labels=2,
        below_label=" ",
        above_label=" ",
        vertical=False,
        width=0.3,
        height=0.05,
        position_x=0.025,
    )
    if add_mesh_kws is not None and "scalar_bar_args" in add_mesh_kws:
        default_scalar_bar_args = {**default_scalar_bar_args, **add_mesh_kws["scalar_bar_args"]}

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
        "name": "mesh",
        "opacity": 1.0,
    }

    # combine the default and user-given kwargs
    default_add_mesh_kws = default_add_mesh_kws if add_mesh_kws is None else {**default_add_mesh_kws, **add_mesh_kws}
    default_add_mesh_kws["scalar_bar_args"] = default_scalar_bar_args

    plotter.add_mesh(
        mesh=mesh,
        scalars=field,
        **default_add_mesh_kws,
    )

    if free_boundaries:
        draw_free_boundaries(
            get_free_boundaries(mesh),
            plotter=plotter
        )

    return plotter


def plot_electrograms(
    times,
    electrograms,
    names=None,
    woi=None,
    y_separation=1,
    y_start=0,
    colour=None,
    axes=None,
):
    """
    Plot electrogram traces.

    Args:
        times (ndarray): times at which voltages were measured
        electrograms (ndarray): Electrogram traces. Two-dimensional of size N_points x N_times for bipolar voltage,
            or two-dimensional of shape N_points x N_times x 2 for unipolar dimensional.
        names (ndarray, optional): List of electrode names, on per electrogram. If provided, names these will be used to
            label each electrogram.
        woi (tuple, optional): start and stop times of the window of interest. If provided, dashed vertical lines
            will be added to the plot at these times.
        y_separation (float, optional): Vertical spacing to add between consecutive electrograms.
        y_start (float, optional): The first electrogram will have this value added to it (to shift the electrogram
            up or down the y axis).
        colour (str or list, optional): Colour or list of colours to use for plotting.
        axis (matplotlib.axes.Axes, optional): Matplotlib Axes on which to plot the traces. If None, a new figure and axes
            will be created.

    Returns:
        figure (matplotlib.Figure): Figure on which the traces have been plotted
        axis (matplotlib.axes.Axes): Axes on which the traces have been plotted.
    """

    separations = y_start + np.arange(electrograms.shape[0]) * y_separation
    colour = "xkcd:cerulean" if colour is None else colour

    if axes is None:
        figure, axes = plt.subplots(constrained_layout=True, figsize=(6, 0.4*len(electrograms)))
    else:
        figure = axes.get_figure()

    # Plot electrograms
    axes.plot(times, electrograms.T + separations, label=names, color=colour)

    # Add names
    if names is not None:
        axes.set_yticks(separations)
        axes.set_yticklabels(names)

    # Add a horizontal line for each electrogram at its zero voltage position
    for y in separations:
        axes.axhline(y, color='grey', linestyle='--', linewidth=0.8, alpha=0.6)

    # Vertical line at the window of interest
    if woi is not None:
        woi_start, woi_stop = woi
        axes.axvline(woi_start, color="grey", linestyle='--', linewidth=0.8, alpha=0.6)
        axes.axvline(woi_stop, color="grey", linestyle='--', linewidth=0.8, alpha=0.6)

    # Remove the border and ticks
    plt.tick_params(axis='both', which='both', length=0)
    for spine in ['left', 'right', 'top']:
        axes.spines[spine].set_visible(False)

    return figure, axes
