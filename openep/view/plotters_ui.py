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
        
from PyQt5 import  QtWidgets, QtGui


def create_plotter_layout(plotter):
    """Create the layout of a BackgroundPlotter.

    Args:
        plotter [BackgrounPlotter]: Plotter for which the layout will be created

    Returns:
        plotter_layout [QtWidgets.QVBoxLayout]: Layout for the plotter
    """

    # TODO: We're currently using line edits to set colourbar limits - look into RangeSlider
    plotter_layout = QtWidgets.QVBoxLayout(plotter)
    plotter.colour_bar_layout, plotter.lower_limit, plotter.upper_limit = _create_colourbar_layout()
    plotter.opacity_layout, plotter.opacity = _create_opacity_layout()
    control_layout = _create_control_layout(plotter=plotter)
    plotter_layout.addLayout(control_layout)
    plotter_layout.addStretch()

    return plotter_layout

def _create_colourbar_layout():
    """Create a layout with widgets for setting the limits of the colourbar for a BackgroundPlotter."""

    colour_bar_layout = QtWidgets.QHBoxLayout()

    colour_bar_text = QtWidgets.QLabel("Colourbar limits:")
    colour_bar_text.setMinimumWidth(100)
    colour_bar_text.setMaximumWidth(200)
    colour_bar_text.setStyleSheet('QLabel {border: 0px; background-color: #d8dcd6;}')
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
    opacity_layout.addStretch

    opacity_text = QtWidgets.QLabel("Opacity:")
    opacity_text.setMinimumWidth(50)
    opacity_text.setMaximumWidth(120)
    opacity_text.setStyleSheet('QLabel {border: 0px; background-color: #d8dcd6;}')
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

def _create_control_layout(plotter):
    """Add the colourbar and opacity layouts to a single horizontal layout"""

    control_layout = QtWidgets.QHBoxLayout()
    control_layout.addLayout(plotter.colour_bar_layout)
    control_layout.addLayout(plotter.opacity_layout)

    return control_layout
