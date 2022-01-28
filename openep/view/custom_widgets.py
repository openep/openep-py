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

from PySide2 import QtCore, QtWidgets, QtGui
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar


class CustomDockWidget(QtWidgets.QDockWidget):
    """
    A draggable and resizable dockwidget.

    From: https://stackoverflow.com/a/63642583
    """

    def __init__(self, title: str):

        super().__init__(title)
        self._init_title()
        self.dockLocationChanged.connect(self.on_dockLocationChanged)

    def _init_title(self):

        self.setTitleBarWidget(QtWidgets.QWidget())
        title_font = QtGui.QFont()
        title_font.setBold(True)
        self.setFont(title_font)

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
    """A pared-down matplotlib toolbar."""

    def __init__(self, canvas_, parent_, keep_actions=None):

        super().__init__(canvas_, parent_)

        self._keep_actions = keep_actions if keep_actions is not None else ['Home', 'Back', 'Forward', 'Zoom', 'Pan', 'Save']
        self._remove_unwanted_buttons()

    def _remove_unwanted_buttons(self):

        # the following snippet is from: https://stackoverflow.com/a/63341907
        for action in self.actions()[:-1]:  # we always want to display the xy coordinate
            if action.text() not in self._keep_actions:
                self.removeAction(action)
