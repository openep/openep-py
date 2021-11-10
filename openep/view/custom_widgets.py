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
Custom Qt Widgets for the OpenEP-Py GUI.
"""

from PyQt5 import QtCore, QtWidgets
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar


class CustomDockWidget(QtWidgets.QDockWidget):
    def __init__(self, title: str):
        super().__init__(title)
        self.setTitleBarWidget(QtWidgets.QWidget())
        self.dockLocationChanged.connect(self.on_dockLocationChanged)

    def on_dockLocationChanged(self):
        main: QtWidgets.QMainWindow = self.parent()
        all_dock_widgets = main.findChildren(QtWidgets.QDockWidget)

        for dock_widget in all_dock_widgets:
            sibling_tabs = main.tabifiedDockWidgets(dock_widget)
            # If you pull a tab out of a group the other tabs still see it as a sibling while dragging...
            sibling_tabs = [s for s in sibling_tabs if not s.isFloating()]

            if len(sibling_tabs) != 0:
                # Hide title bar
                dock_widget.setTitleBarWidget(QtWidgets.QWidget())
            else:
                # Re-enable title bar
                dock_widget.setTitleBarWidget(None)

    def minimumSizeHint(self) -> QtCore.QSize:
        return QtCore.QSize(100, 100)


class CustomNavigationToolbar(NavigationToolbar):

    def __init__(self, canvas_, parent_):
        
        super().__init__(canvas_, parent_)
        self._remove_unwanted_actions()

    def _remove_unwanted_actions(self):

        # the following snippet is from: https://stackoverflow.com/a/63341907
        unwanted_buttons = ['Pan', 'Subplots', 'Customize']
        for action in self.actions():
            if action.text() in unwanted_buttons:
                self.removeAction(action)
