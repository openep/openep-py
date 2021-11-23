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
Create and manipulate Matplotlib canvases.
"""

from PyQt5 import QtWidgets
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.widgets
import matplotlib.pyplot as plt

import openep.view.custom_widgets


def create_canvas():
    """
    Create an interactive matploblib canvas.
    """

    figure, axis = plt.subplots(ncols=1, nrows=1)
    figure.set_facecolor("white")

    # hide it until we have data to plot
    axis.axis('off')
    # and don't display xy coordinates in the toolbar when hovering over the axis
    axis.format_coord = lambda x, y: ""

    canvas = FigureCanvas(figure)

    return canvas, figure, axis


def create_egm_selection_layout():
    """Create a layout with widgets for selecting which electrograms to plot."""

    egm_selection_layout = QtWidgets.QHBoxLayout()

    egm_selection_text = QtWidgets.QLabel("Select EGMs (indices of points):")
    egm_selection_text.setMinimumWidth(220)
    egm_selection_text.setMaximumWidth(300)
    egm_selection_text.setStyleSheet('border: 0px; background-color: white;')

    egm_selection = QtWidgets.QLineEdit()
    egm_selection.setMinimumWidth(100)
    egm_selection.setMaximumWidth(300)
    egm_selection.setStyleSheet('border: 1px solid #d8dcd6; background-color: white;')
    egm_selection.setText("0")
    egm_selection.setPlaceholderText("indices")

    egm_selection_layout.addWidget(egm_selection_text)
    egm_selection_layout.addWidget(egm_selection)
    egm_selection_layout.addStretch()

    return egm_selection_layout, egm_selection


def create_egm_type_layout():
    """Create a layout with widgets for selecting which type(s) of electrograms to plot."""

    reference_checkbox = QtWidgets.QCheckBox("Reference")
    reference_checkbox.setStyleSheet("color: #be0119; background-color: white;")  # xkcd:scarlet
    reference_checkbox.setChecked(False)

    bipolar_checkbox = QtWidgets.QCheckBox("Bipolar")
    bipolar_checkbox.setStyleSheet("color: #0485d1; background-color: white;")  # xkcd:cerulean
    bipolar_checkbox.setChecked(True)

    unipolar_A_checkbox = QtWidgets.QCheckBox("Unipolar: A")
    unipolar_A_checkbox.setStyleSheet("color: #2a7e19; background-color: white;")  # xkcd:tree green
    unipolar_A_checkbox.setChecked(True)

    unipolar_B_checkbox = QtWidgets.QCheckBox("Unipolar: B")
    unipolar_B_checkbox.setStyleSheet("color: #fb7d07; background-color: white;")  # xkcd:pumpkin
    unipolar_B_checkbox.setChecked(True)

    egm_type_layout = QtWidgets.QHBoxLayout()
    egm_type_layout.addWidget(reference_checkbox)
    egm_type_layout.addWidget(bipolar_checkbox)
    egm_type_layout.addWidget(unipolar_A_checkbox)
    egm_type_layout.addWidget(unipolar_B_checkbox)
    egm_type_layout.addStretch()

    return egm_type_layout, reference_checkbox, bipolar_checkbox, unipolar_A_checkbox, unipolar_B_checkbox


def add_toolbar(canvas, parent, keep_actions=None):
    """
    Add a toolbar to a canvas.

    Args:
        keep_actions (list): list of actions to be kept in the toolbar.
            The default is None, in which case the following actions
            will be kept: 'Home', 'Back', 'Forward', 'Zoom', 'Save'.
    """

    toolbar = openep.view.custom_widgets.CustomNavigationToolbar(
        canvas,
        parent,
        keep_actions=keep_actions,
    )

    return toolbar


def create_canvas_widget(canvas, toolbar):
    """Create a placeholder widget to hold a toolbar and canvas.
    """

    canvas_layout = QtWidgets.QVBoxLayout()
    canvas_layout.addWidget(canvas)
    canvas_layout.addWidget(toolbar)

    canvas_widget = QtWidgets.QWidget()
    canvas_widget.setLayout(canvas_layout)
    canvas_widget.setStyleSheet("border-width: 0px; border: 0px; background-color: #d8dcd6;")

    return canvas_widget


def add_woi_range_slider(figure, axis, valmin, valmax, valstep=5):
    """
    Add a RangeSlider to a figure. The RangeSlider is for setting the window of interest.

    Args:
        figure (matplotlib.Figure): figure to which the slider will be added.
        axis (list): list of points ([xstart, ystart, xlength, ylength]) that will be
            used to define the position of the RangeSlider on the figure.
        valmin ([type]): Minimum value of the slider.
        valmax ([type]): Maximum value of the slider.
        valstep (int, optional): The slider can take values in this step size. Defaults to 5.
    """

    slider_axis = figure.add_axes(axis)
    slider = matplotlib.widgets.RangeSlider(
        ax=slider_axis,
        label="WOI",
        valmin=valmin,
        valmax=valmax,
        closedmin=True,
        closedmax=True,
        dragging=True,
        valstep=valstep,
        facecolor="xkcd:light grey",
    )

    return slider


def add_woi_set_button(figure, axis):
    """
    Add a button to set the window of interest.

    When pressed, the electrogram traces will be used to interpolate
    data onto the surface of the mesh and draw the map if nececssary.

    Args:
        figure (matplotlib.Figure): figure to which the button will be added.
        axis (list): list of points ([xstart, ystart, xlength, ylength]) that will be
            used to define the position of the Button on the figure.
    """

    set_woi_axis = figure.add_axes(axis)
    set_woi_button = matplotlib.widgets.Button(
        ax=set_woi_axis,
        label="Set WOI",
        color="xkcd:light grey",
    )

    return set_woi_button
