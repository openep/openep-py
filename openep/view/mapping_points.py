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
Class and functions for handling and displaying list of mapping points.
"""

from typing import Optional

from PySide2 import QtCore, QtWidgets, QtGui
from PySide2.QtCore import Qt

import numpy as np

from openep.view.custom_widgets import CustomDockWidget
from openep.view.system_manager import System

class MappingPointsModel(QtCore.QAbstractTableModel):
    """Table model for handling data associated with mapping points.
    
    The following data for each mapping point is stored in a 2D Numpy Array:
        - index (0-based)
        - tag  (physician-visible name of the point applied during the clinical case)
        - name (internal name used by clinical mapping system)
        - bipolar voltage
        - local activation time

    """

    def __init__(
        self,
        system,
    ):
            
        super().__init__()
        self._system = system
        self._init_headers()
        self._init_data()
    
    @property
    def system(self):
        return self._system
    
    @system.setter
    def system(self, system):
        self._system = system
        self._init_data()

    def _init_headers(self):
        
        self.headers = ["Index", "Tag", "Name", "Voltage", "LAT"]
        
        # these will later be used for showing/hiding specific columns
        self._column_indices = {header: index for index, header in enumerate(self.headers)}

    def headerData(self, section, orientation, role):
        """Set the table header."""

        if role == QtCore.Qt.DisplayRole :
            if orientation == QtCore.Qt.Horizontal :
                return self.headers[section]
            elif orientation == QtCore.Qt.Vertical :
                return section
        return None

    def _init_data(self):
        """Extract mapping point data from the system.
        
        Store data in a 2d array, as required for displaying in a table.
        """
        
        if self._system is None:
            self._data = None
            self.include = np.array([])
            return

        electric = self._system.case.electric
        has_bipolar = True if electric.bipolar_egm.egm is not None else False
        
        if not has_bipolar:
            self._data = None
            self.include = np.array([])
            return

        self.include = electric.include
        self._index = np.arange(electric.bipolar_egm.voltage.shape[0], dtype=int)
        
        # TODO: should the activation time be absolute, or relative to the window of interest?
        activation_time = electric.annotations.local_activation_time
        self._data = np.concatenate(
            [
                self._index[:, np.newaxis],
                electric.names[:, np.newaxis],
                electric.internal_names[:, np.newaxis],
                electric.bipolar_egm.voltage[:, np.newaxis],
                activation_time[:, np.newaxis],
            ],
            axis=1,
        )
        self.layoutChanged.emit()

    def data(self, index, role):
        
        if role == Qt.DisplayRole:

            value = self._data[index.row(), index.column()]
            return value

    def rowCount(self, index=None):
        return 0 if self._data is None else self._data.shape[0]

    def columnCount(self, index=None):
        return 0 if self._data is None else self._data.shape[1]


class SortProxyModel(QtCore.QSortFilterProxyModel):
    """Model for sorting and filtering data in a table."""
    
    def lessThan(self, left_index, right_index):
        """Compare two neighbouring values, for sorting rows."""

        left_var = self.sourceModel().data(left_index, Qt.DisplayRole)
        right_var = self.sourceModel().data(right_index, Qt.DisplayRole)
        
        try:
            return float(left_var) < float(right_var)
        except ValueError as e:
            if not str(e).startswith("could not convert string to float:"):
                raise e
            return left_var < right_var

class MappingPointsDock(CustomDockWidget):
    """A dockable widget for handling and displaying the mapping points data."""

    def __init__(self, title: str, system: Optional[System] = None, model=None, proxy_model=None):

        super().__init__(title)
        
        self.model = model if model is not None else MappingPointsModel(system=system)
        self.proxy_model = proxy_model if proxy_model is not None else SortProxyModel()
        self.proxy_model.setSourceModel(self.model)

        self.table = QtWidgets.QTableView()
        self.table.setModel(self.proxy_model)
        self.table.verticalHeader().hide()
        self.table.setSortingEnabled(True)
        self.table.resizeColumnsToContents()
        self.table.setShowGrid(False)
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.table.setAlternatingRowColors(True)
        self.table.setStyleSheet(
            "alternate-background-color: #262E38;"
        )
        
        self.main = QtWidgets.QMainWindow()
        self.main.setCentralWidget(self.table)
        self.setWidget(self.main)
        
        self._init_header_context_menu()
        self._init_table_context_menu()

    def _init_header_context_menu(self):
        """Create menu for selecting which columns to show when right-clicking on the header"""

        header = self.table.horizontalHeader()
        header.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        header.customContextMenuRequested.connect(self._launch_header_context_menu)
        
        self.header_menu = QtWidgets.QMenu()
        show_menu = QtWidgets.QMenu("Show/hide")
        self.header_menu.addMenu(show_menu)
        
        for column_name in self.model.headers:
            action = QtWidgets.QAction(column_name, self.main, checkable=True, checked=True)
            show_menu.addAction(action)
            action.toggled.connect(
                lambda checked, name=column_name: self.column_visibility(checked, self.model._column_indices[name])
            )
        
    def _launch_header_context_menu(self, pos):
        self.header_menu.exec_(QtGui.QCursor.pos())

    def column_visibility(self, show, column_index):
        self.table.setColumnHidden(column_index, not show)

    def _init_table_context_menu(self):
        """Create a menu for hiding selected rows."""

        self.table.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.table.customContextMenuRequested.connect(self._launch_table_context_menu)

        self.table_menu = QtWidgets.QMenu()
        self.delete_points = QtWidgets.QAction("Delete point(s)", self.main, checkable=False)
        self.table_menu.addAction(self.delete_points)
        
    def _launch_table_context_menu(self, pos):
        self.table_menu.exec_(QtGui.QCursor.pos())


class RecycleBinDock(MappingPointsDock):
    """A dockable widget for handling and displaying deleted mapping points."""

    def __init__(self, title: str, system: Optional[System] = None, model=None, proxy_model=None):

        super().__init__(title, system, model, proxy_model)

    def _init_table_context_menu(self):
        """Create a menu for hiding selected rows. These """

        self.table.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.table.customContextMenuRequested.connect(self._launch_table_context_menu)

        self.table_menu = QtWidgets.QMenu()
        self.restore_points = QtWidgets.QAction("Restore point(s)", self.main, checkable=False)
        self.table_menu.addAction(self.restore_points)
