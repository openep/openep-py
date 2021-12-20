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
Functions for creating widgets for individual systems
"""

from PyQt5.QtCore import Qt
from PyQt5 import QtGui, QtWidgets
from .custom_widgets import CustomDockWidget

def create_system_dock(plotter):
    """Create a new dockable widget that will contain a pyvista-qt plotter for rendering 3D maps.

    Args:
        plotter (BackgroundPlotter): plotter to set at the widget for the new dock

    Returns:
        dock (CustomDockWidget): dockable widget containing the plotter
    """

    dock = CustomDockWidget("Temporary title")  # this will be changed
    dock.main = QtWidgets.QMainWindow()

    # The dock is set to have bold font (so the title stands out)
    # But all other widgets should have normal weighted font
    main_font = QtGui.QFont()
    main_font.setBold(False)
    dock.main.setFont(main_font)

    # Create a menubar for this dock
    # From here we can control e.g. the scalar values projected onto the surface
    dock.main.menubar = dock.main.menuBar()
    dock.main.menubar.setNativeMenuBar(False)
    dock.main.menubar

    # Set widget
    dock.main.setCentralWidget(plotter)
    dock.setWidget(dock.main)
    dock.setAllowedAreas(Qt.AllDockWidgetAreas)

    return dock

def add_field_menu(dock, plotter, system_name, scalar_fields):
    """Add a Field menu to the menubar. This is used for selecting the scalar field to project onto the surface.

    Args:
        dock (CustomDockWidget): the dockwidget to which the 'Field' menu will be added to the menubar
        plotter (BackgroundPlotter): The plotter which will have actions added for when specific values from the Field
            menu are selected.
        system_name (str): Name of the system. Will be used in setting the title of the dock widget.
        scalar_fields (dict): Names and values of the scalar fields to be added as options in the Field menu.
    """

    field_menu = dock.main.menubar.addMenu("Field")
    field_group = QtWidgets.QActionGroup(dock.main)

    plotter.scalar_field_actions = {}
    for field_name in scalar_fields:

        dock.setWindowTitle(f"{system_name}: {field_name}")
        action = QtWidgets.QAction(field_name, dock.main, checkable=True)
        action.setChecked(False)

        field_menu.addAction(action)
        field_group.addAction(action)
        plotter.scalar_field_actions[field_name] = action

    return dock, plotter
