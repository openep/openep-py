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
Class for creating the system manager widget.
"""
from PySide2 import QtGui, QtWidgets
from PySide2.QtWidgets import QAction
from PySide2.QtCore import Qt
from .custom_widgets import CustomDockWidget


class SystemManagerDockWidget(CustomDockWidget):
    """A dockable widget that for managing systems loaded into the GUI"""

    def __init__(self, title: str):

        super().__init__(title)
        self._init_main_window()
        self._create_system_manager_menubar()
        self._create_system_manager_table()

        # Setting nested layouts
        self.main.setCentralWidget(self.table)
        self.setWidget(self.main)

    def _init_main_window(self):
        """Create a MainWindow that will be the widget of the system manager"""

        self.main = QtWidgets.QMainWindow()

    def _create_system_manager_menubar(self):
        """Add a menu bar to the System Manager.

        The file menu has options for loading either an OpenEP dataset or an openCARP one, as well as
        loading auxillary data into an existing openCARP system.

        The view menu has options for creating a new 3D-viewer for an existing system.
        """

        # For loading new systems, deleting systems, adding a view of existing systems
        # Adding data to systems, exporting system to a different format
        self.main.menubar = self.main.menuBar()
        self.main.menubar.setNativeMenuBar(False)

        file_menu = self.main.menubar.addMenu("File")

        load_menu = QtWidgets.QMenu("Load", self.main)
        load_openep_mat_action = QAction("OpenEP", self.main)
        load_menu.addAction(load_openep_mat_action)
        load_opencarp_action = QAction("openCARP", self.main)
        load_menu.addAction(load_opencarp_action)

        self.main.load_openep_mat_action = load_openep_mat_action
        self.main.load_opencarp_action = load_opencarp_action

        # The below will be used for adding data to openCARP datasets
        self.main.add_data_menu = QtWidgets.QMenu("Add data to", self.main)

        # The below will be used exporting OpenEP mesh data in openCARP format
        self.main.export_data_menu = QtWidgets.QMenu("Export", self.main)

        file_menu.addMenu(load_menu)
        file_menu.addMenu(self.main.add_data_menu)
        file_menu.addSeparator()
        file_menu.addMenu(self.main.export_data_menu)

        # The below will be used for creating a new 3D viewer for a given system
        view_menu = self.main.menubar.addMenu("View")
        self.main.add_view_menu = QtWidgets.QMenu("Add view for", self.main)
        view_menu.addMenu(self.main.add_view_menu)
        view_menu.addSeparator()

    def _create_system_manager_table(self):
        """
        Create a widget that will display information about each system loaded, and
        allow selection of an active system.
        """

        self.table = QtWidgets.QWidget()
        self.table_layout = QtWidgets.QGridLayout()
        self.table.setLayout(self.table_layout)
        
        # We're going to add an additional column to take care of spacing
        self.table_layout.setHorizontalSpacing(0)
        
        # We use a large number here so the rows do not change size to fill the vertical space
        self.table_layout.setRowStretch(1e6, 1)

        # add a row containing heading labels to the grid
        self._add_heading_bar()

        # we need a method for selecting the active system
        self._active_system_button_group = QtWidgets.QButtonGroup()

    def _add_heading_bar(self):
        """Add a heading bar to the table layout."""

        heading_font = QtGui.QFont()
        heading_font.setBold(True)

        heading_names = [
            'Name',
            'File',
            'Type',
            'Active',
        ]
        self._alignments = [
            Qt.AlignLeft,
            Qt.AlignLeft,
            Qt.AlignLeft,
            Qt.AlignCenter,
        ]
        self._column_stretches = [3, 8, 2, 1]
        row = 0
        row_width = 1
        
        for column, heading_name in enumerate(heading_names):
            
            heading = QtWidgets.QLabel(heading_name)
            heading.setStyleSheet('background: #60798B; padding: 5px;')
            heading.setFont(heading_font)

            # For the first 3 columns, add an additional column to the right hand side.
            # This adds a small space between the columns.
            if column < len(heading_names) - 1:
                self.table_layout.addWidget(heading, row, 2 * column, row_width, column_width:=2)
                self.table_layout.setColumnMinimumWidth(2 * column + 1, 20 * self._column_stretches[column])
                self.table_layout.setColumnStretch(2 * column, self._column_stretches[column])
            else:
                self.table_layout.addWidget(heading, row, 2 * column, row_width, column_width:=1, self._alignments[column])

    def create_export_action(self, system_basename, export_name):
        """Add an action to the Menubar for exporting a system dataset to a specific format.

        Args:
            system_basename (str): Basename of the system to be exported (i.e. system.basename)
            export_name (str): Name of the export action action (will appear under 'File > export > `system.basename`')

        Returns:
            export_action (QAction): Action to perform for exporting the data
        """

        # TODO: the menu should be created separately from the action
        export_menu = QtWidgets.QMenu(system_basename, self)
        export_action = QAction(export_name, self)
        export_menu.addAction(export_action)

        self.main.export_data_menu.addMenu(export_menu)

        return export_action

    def create_view_action(self, system_basename):
        """Add an action to the Menubar for creating a new 3d viewer of a system.

        Args:
            system_basename (str): Basename of the system to be exported (i.e. system.basename)

        Returns:
            view_action (QAction): Action to open a new 3d viewer
        """

        view_action = QAction(system_basename, self)
        self.main.add_view_menu.addAction(view_action)

        return view_action

    def create_add_data_action(self, system_basename, data_type):
        """Add an action to the menubar for loading auxillary data into a system.

        Args:
            system_basename (str): Basename of the system to be exported (i.e. system.basename)
            data_type (str): Description of data that will be exported using this action
        """

        # TODO: the menu should be created separately from the action
        add_data_menu = QtWidgets.QMenu(system_basename, self)
        add_data_action = QAction(data_type, self)
        add_data_menu.addAction(add_data_action)
        self.main.add_data_menu.addMenu(add_data_menu)

        return add_data_action

    def update_system_manager_table(
        self,
        name: str,
        basename: str,
        data_type: str,
        is_active: str,
        row_number: int,
    ):
        """Update the System manager table when a new system is loaded.

        We need to remove the vertical spacer, add a new row in the system table,
        then put the vertical spacer back.

        We also need to add the system to the 'File > Add data to' submenu and the
        'View > Add view for' submenu.

        Args:
            name (str): label for the system
            basename (str): basename of the file(s) for the given system
            data_type (str): type of system data - either OpenEP or openCARP
            is_active (bool): whether the system is currently the active system in the GUI
            row_number (int): which row number in the QGridLayout to add the new row

        Returns:
            name_widget (QtWidgets.QLineEdit): widget for editing the name of the system
            basename_widget (QtWidgets.QLabel): widget displaying the basename of the system
            data_type_widget (QtWidgets.QLabel): widget displaying the data type of the system
            active_widget (QtWidgets.QRadioButton): widget for setting the system to be the active system
        """

        widgets = self._create_widgets_for_system_row(
            name=name,
            basename=basename,
            data_type=data_type,
            is_active=is_active,
            )
        name_widget, basename_widget, data_type_widget, active_widget = widgets
        self._active_system_button_group.addButton(active_widget)
        
        row_width = 1        
        for column, widget in enumerate(widgets):

            # For the first 3 columns, add an additional column to the right hand side.
            # This adds a small space between the columns.
            if column < len(widgets) - 1:
                self.table_layout.addWidget(widget, row_number, 2 * column, row_width, column_width:=2)
                self.table_layout.setColumnMinimumWidth(2 * column + 1, 20 * self._column_stretches[column])
                self.table_layout.setColumnStretch(2 * column, self._column_stretches[column])
            else:
                self.table_layout.addWidget(widget, row_number, 2 * column, row_width, column_width:=1, self._alignments[column])

        return name_widget, basename_widget, data_type_widget, active_widget

    def _create_widgets_for_system_row(
        self,
        name: str,
        basename: str,
        data_type: str,
        is_active: str,
    ):
        """Create widgets with info about a given system.

        Args:
            name (str): label for the system
            basename (str): basename of the file(s) for the given system
            data_type (str): type of system data - either OpenEP or openCARP
            is_active (bool): whether the system is currently the active system in the GUI

        Returns:
            name_widget (QtWidgets.QLineEdit): widget for editing the name of the system
            basename_widget (QtWidgets.QLabel): widget displaying the basename of the system
            data_type_widget (QtWidgets.QLabel): widget displaying the data type of the system
            active_widget (QtWidgets.QRadioButton): widget for setting the system to be the active system
        """
       
        name_widget = QtWidgets.QLineEdit()
        name_widget.setText(str(name))
        name_widget.setPlaceholderText("system name")
        
        basename_widget = QtWidgets.QLabel()
        basename_widget.setText(str(basename))
        
        data_type_widget = QtWidgets.QLabel()
        data_type_widget.setText(data_type)
        
        active_widget = QtWidgets.QRadioButton()
        active_widget.setChecked(is_active)

        return name_widget, basename_widget, data_type_widget, active_widget
