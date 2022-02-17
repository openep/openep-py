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

    def headerData(self, section, orientation, role):
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
            return

        electric = self._system.case.electric
        has_bipolar = True if electric.bipolar_egm.egm is not None else False
        
        if not has_bipolar:
            self._data = None
            return

        if not hasattr(electric, '_include'):
            electric._include = np.full_like(
                electric.bipolar_egm.voltage,
                fill_value=True,
                dtype=bool,
            )

        activation_time = electric.annotations.local_activation_time - electric.annotations.reference_activation_time
        self._data = np.concatenate(
            [
                np.arange(electric.bipolar_egm.voltage.shape[0], dtype=int)[:, np.newaxis],
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
    
    def lessThan(self, left_index, right_index):

        left_var = self.sourceModel().data(left_index, Qt.DisplayRole)
        right_var = self.sourceModel().data(right_index, Qt.DisplayRole)
        
        try:
            return float(left_var) < float(right_var)
        except ValueError as e:
            if not str(e).startswith("could not convert string to float:"):
                raise e
            return left_var < right_var

class MappingPointsDock(CustomDockWidget):
    """A dockable widget for handling and viewing the mapping points data."""

    def __init__(self, title: str, system: Optional[System] = None):

        super().__init__(title)
        
        self.model = MappingPointsModel(system=system)
        self.proxy_model = SortProxyModel()
        self.proxy_model.setSourceModel(self.model)
        self.table = QtWidgets.QTableView()
        self.table.setModel(self.proxy_model)
        self.table.verticalHeader().hide()
        self.table.setSortingEnabled(True)
        self.table.setShowGrid(False)
        self.table.setAlternatingRowColors(True)
        self.table.setStyleSheet(
            "alternate-background-color: #262E38;"
        )

        header = self.table.horizontalHeader()
        header.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        header.customContextMenuRequested.connect(self.header_context_menu)
        
        self.main = QtWidgets.QMainWindow()
        self.main.setCentralWidget(self.table)
        self.setWidget(self.main)

    def header_context_menu(self, pos):
        
        # column_index = self.table.horizontalHeader().logicalIndexAt(pos)
        
        menu = QtWidgets.QMenu()
  
        show_menu = QtWidgets.QMenu("Show/hide")
        show_index = QtWidgets.QAction("Index")
        show_tag = QtWidgets.QAction("Tag")
        show_name = QtWidgets.QAction("Name")
        show_voltage = QtWidgets.QAction("Voltage")
        show_lat = QtWidgets.QAction("LAT")
        show_menu.addActions(
            [
                show_index,
                show_tag,
                show_name,
                show_voltage,
                show_lat,
            ]
        )
        menu.addMenu(show_menu)

        menu.exec_(QtGui.QCursor.pos())
