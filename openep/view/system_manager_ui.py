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
Functions for creating the system manager widget.
"""
from typing import Union
from attr import attrs

from PyQt5 import QtGui, QtWidgets
from PyQt5.QtCore import Qt

import openep
import openep.view.custom_widgets
import openep.view.plotters


def add_heading_bar(layout):
    """Add a heading bar to the system managet layout.

    Args:
        layout (QGridLayout): Layout to which the heading bar will be added.
    """

    heading_font = QtGui.QFont()
    heading_font.setBold(True)

    heading_bar = QtWidgets.QWidget()
    heading_bar_layout = QtWidgets.QHBoxLayout()
    total_width = 0
    for heading_name, width in zip(
        ['Name', 'Directory', 'Type', 'Active'],
        [2, 8, 1, 1],
    ):
        heading = QtWidgets.QLabel(heading_name)
        heading.setFont(heading_font)
        heading.setFrameShape(QtWidgets.QFrame.NoFrame)
        heading.setLineWidth(0)
        heading_bar_layout.addWidget(heading, width)
        total_width += width

    heading_bar.setLayout(heading_bar_layout)
    heading_bar.setStyleSheet('QWidget {background-color: #95a3a6;}')
    layout.addWidget(heading_bar, 0, 0, total_width, 1)


def create_vertical_stretch():
    """Create a vertical spacer widget"""

    vertical_stretch = QtWidgets.QSpacerItem(
        0, 0,
        QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding,
    )

    return vertical_stretch


def create_widgets_for_system_row(
    layout: QtWidgets.QGridLayout,
    name: str,
    directory: str,
    data_type: str,
    is_active: str,
    row_number: int,
    ):
    """Create a QHBoxLayout containing widgets with info about a given system.

    The row is added to the layout in place. The created widgets are also returned.

    Args:
        layout (QGridLayout): Layout to which the heading bar will be added.
        name (str): label for the system
        directory (str): path to the directory contianing the files for the given system
        data_type (str): type of system data - either OpenEP or openCARP
        is_active (bool): whether the system is currently the active system in the GUI
        row_number (int): which row number in the QGridLayout to add the new row

    Returns:
        name_widget (QtWidgets.QLineEdit): widget for editing the name of the system
        directory_widget (QtWidgets.QLabel): widget displaying the directory of the system
        data_type_widget (QtWidgets.QLabel): widget displaying the data type of the system
        active_widget (QtWidgets.QRadioButton): widget for setting the system to be the active system
    """

    system_row = QtWidgets.QWidget()
    system_layout = QtWidgets.QHBoxLayout()

    # The widths of the widgets are the same as those used in add_heading_bar()
    name_width, directory_width, data_type_width, active_width = [2, 8, 1, 1]
    total_width = 12

    name_widget = QtWidgets.QLineEdit()
    name_widget.setText(str(name))
    name_widget.setPlaceholderText("system name")
    system_layout.addWidget(name_widget, name_width)

    directory_widget = QtWidgets.QLabel()
    directory_widget.setText(str(directory))
    system_layout.addWidget(directory_widget, directory_width)

    data_type_widget = QtWidgets.QLabel()
    data_type_widget.setText(data_type)
    system_layout.addWidget(data_type_widget, data_type_width)

    active_widget = QtWidgets.QRadioButton()
    active_widget.setChecked(is_active)
    system_layout.addWidget(active_widget, active_width)

    system_row.setLayout(system_layout)
    system_row.setStyleSheet('background-color: white;')

    layout.addWidget(system_row, row_number, 0, total_width, 1)

    return name_widget, directory_widget, data_type_widget, active_widget
