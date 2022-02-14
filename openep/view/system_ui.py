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

from PySide2.QtCore import Qt
from PySide2 import QtGui, QtWidgets
from PySide2.QtWidgets import QAction, QActionGroup
from .custom_widgets import CustomDockWidget


def create_system_dock(central_widget):
    """Create a new dockable widget that will contain a pyvista-qt plotter for rendering 3D maps.

    Args:
        central_widget (QWidget): Widget containing the BackgroundPlotter to add to the new dock

    Returns:
        dock (CustomDockWidget): dockable widget containing the plotter
    """

    dock = CustomDockWidget("Temporary title")  # this will be changed
    dock.main = QtWidgets.QMainWindow()

    # Create a menubar for this dock
    # From here we can control e.g. the scalar values projected onto the surface
    dock.main.menubar = dock.main.menuBar()
    dock.main.menubar.setNativeMenuBar(False)
    dock.main.menubar

    # Set widget
    dock.main.setCentralWidget(central_widget)
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

    Returns:
        dock (CustomDockWidget): Dockwidget with a 'Field' menu added to the menubar.
        plotter (BackgroundPlotter): Plotter with actions for selecting the scalar field stored in a dictionary
            as the attribute `plotter.scalar_field_actions`.
    """

    field_menu = dock.main.menubar.addMenu("Field")
    field_group = QActionGroup(dock.main)

    plotter.scalar_field_actions = {}
    for field_name in scalar_fields:

        dock.setWindowTitle(f"{system_name}: {field_name}")
        action = QAction(field_name, dock.main, checkable=True)
        action.setChecked(False)

        field_menu.addAction(action)
        field_group.addAction(action)
        plotter.scalar_field_actions[field_name] = action

    return dock, plotter


def add_show_menu(dock, plotter):
    """Add a Show menu to the menubar. This is used for showing/hiding the mesh, mapping points, surface-projected mapping points.

    Args:
        dock (CustomDockWidget): the dockwidget to which the 'Show' menu will be added to the menubar
        plotter (BackgroundPlotter): The plotter which will have actions added for when specific values from the Show
            menu are selected.

    Returns:
        dock (CustomDockWidget): Dockwidget with a 'Show' menu added to the menubar.
        plotter (BackgroundPlotter): Plotter with actions for showing/hiding actors stored in a dictionary
            as the attribute `plotter.show_actions`.
    """

    show_menu = dock.main.menubar.addMenu("Show/Hide")

    plotter.show_actions = {}
    action_names = ["Surface", "Mapping points", "Surface-projected mapping points"]
    for action_name in action_names:

        action = QAction(action_name, dock.main, checkable=True)
        show = True if action_name == "Surface" else False  # by default only display the mesh
        action.setChecked(show)
        show_menu.addAction(action)
        plotter.show_actions[action_name] = action

    return dock, plotter
