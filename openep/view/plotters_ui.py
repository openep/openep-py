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
Functions for creating the UI for BackgroundPlotters.
"""

from PySide2 import QtCore, QtWidgets, QtGui

import openep.view.static


def create_plotter_layout(plotter, add_link_views_button=True):
    """Create the layout of a BackgroundPlotter.

    Add widgets to the plotter for:
        * setting the colourbar limits (with QLineEdits)
        * setting the opacity of the mesh (with QDoubleSpinBox)
        * linking the view of a secondary plotter with that of the primary plotter (with QPushButton)

    Args:
        plotter [BackgrounPlotter]: Plotter for which the layout will be created

    Returns:
        plotter [BackgrounPlotter]: Plotter for which the layout has been created.
        plotter_layout [QtWidgets.QVBoxLayout]: Layout for the plotter
    """

    # TODO: We're currently using line edits to set colourbar limits - look into RangeSlider
    plotter_layout = QtWidgets.QVBoxLayout(plotter)

    colour_bar_layout, lower_limit, upper_limit = _create_colourbar_layout()
    opacity_layout, opacity = _create_opacity_layout()
    control_layout = _create_control_layout(
        colour_bar_layout=colour_bar_layout,
        opacity_layout=opacity_layout,
    )
    plotter_layout.addLayout(control_layout)

    if add_link_views_button:
        link_views_layout, link_view_with_primary = _create_link_views_layout()
        plotter_layout.addLayout(link_views_layout)
    else:
        link_view_with_primary = None

    plotter_layout.addStretch()

    return plotter_layout, lower_limit, upper_limit, opacity, link_view_with_primary


def _create_colourbar_layout():
    """Create a layout with widgets for setting the limits of the colourbar for a BackgroundPlotter."""

    colour_bar_layout = QtWidgets.QHBoxLayout()

    colour_bar_text = QtWidgets.QLabel("Colourbar limits:")
    colour_bar_text.setMinimumWidth(100)
    colour_bar_text.setMaximumWidth(200)
    #colour_bar_text.setStyleSheet('QLabel {border: 0px; background-color: #d8dcd6;}')
    colour_bar_layout.addWidget(colour_bar_text)

    lower_limit = QtWidgets.QLineEdit()
    lower_limit.setFixedWidth(40)
    lower_limit.setText("0")
    lower_limit.setPlaceholderText("lower")
    lower_limit.setValidator(QtGui.QDoubleValidator(bottom=0.0))
    colour_bar_layout.addWidget(lower_limit)

    upper_limit = QtWidgets.QLineEdit()
    upper_limit.setFixedWidth(40)
    upper_limit.setText("5")
    upper_limit.setPlaceholderText("upper")
    lower_limit.setValidator(QtGui.QDoubleValidator(bottom=0.01))
    colour_bar_layout.addWidget(upper_limit)

    colour_bar_layout.addStretch()

    return colour_bar_layout, lower_limit, upper_limit


def _create_opacity_layout():
    """Create a layout with widgets for setting the opacity of the mesh belonging to a BackgroundPlotter."""

    opacity_layout = QtWidgets.QHBoxLayout()
    opacity_layout.addStretch()

    opacity_text = QtWidgets.QLabel("Opacity:")
    opacity_text.setMinimumWidth(50)
    opacity_text.setMaximumWidth(120)
    #opacity_text.setStyleSheet('QLabel {border: 0px; background-color: #d8dcd6;}')
    opacity_layout.addWidget(opacity_text)

    opacity_selector = QtWidgets.QDoubleSpinBox()
    opacity_selector.setFixedWidth(60)
    opacity_selector.setRange(0.0, 1.0)
    opacity_selector.setSingleStep(0.05)
    opacity_selector.setDecimals(2)
    opacity_selector.setWrapping(False)
    opacity_selector.setValue(1)
    opacity_selector.textFromValue(1)
    opacity_layout.addWidget(opacity_selector)

    return opacity_layout, opacity_selector


def _create_control_layout(colour_bar_layout, opacity_layout):
    """Add colourbar and opacity layouts to a single horizontal layout for controlling mesh properties"""

    control_layout = QtWidgets.QHBoxLayout()
    control_layout.addLayout(colour_bar_layout)
    control_layout.addLayout(opacity_layout)

    return control_layout


def _create_link_views_layout():
    """Add a button to link/unlink the view of a secondary plotter with that of the primary plotter."""

    link_views_layout = QtWidgets.QHBoxLayout()

    link_views_button = QtWidgets.QPushButton()
    link_views_button.setCheckable(True)
    link_views_button.setChecked(True)
    link_views_button.setFixedSize(30, 30)

    link_views_button.setIconSize(QtCore.QSize(30, 30))
    link_views_button.setIcon(QtGui.QIcon(openep.view.static.LINK_ICON))
    link_views_button.setToolTip("3D viewer is linked with the primary 3D viewer")

    link_views_layout.addWidget(link_views_button)
    link_views_layout.addStretch()

    return link_views_layout, link_views_button
