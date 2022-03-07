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
Create a preferences manager.
"""

from PySide2 import QtWidgets, QtCore

__all__ = ["Preferences"]


class PreferencesManager(QtCore.QSettings):
    """Settings manager for OpenEP Gui."""

    widget_mappers = {
        'QCheckBox': ('isChecked', 'setChecked', bool),
        'QLineEdit': ('text', 'setText', str),
        'QSpinBox': ('value', 'setValue', int),
        'QRadioButton': ('isChecked', 'setChecked', bool),
        'QComboBox': ('currentText', 'setCurrentText', str),
        'QSpinBox': ('value', 'setValue', int),
        'QDoublePinBox': ('value', 'setValue', float),
        'QButtonGroup': ('checkedId', 'setId', int),
    }

    settings_changed = QtCore.Signal()

    def __init__(self):
        super().__init__()

        self.settings = QtCore.QSettings("OpenEP-GUI", "settings")
        print(f"Loading preferences from {self.settings.fileName()}")

    def update_widgets_from_settings(self, map):
        for name, widget in map.items():
            cls = widget.__class__.__name__
            getter, setter, dtype = self.widget_mappers.get(cls, (None, None))
            value = self.settings.value(name, type=dtype)
            print("load:", getter, setter, value, type(value), dtype)
            if setter and value is not None:
                fn = getattr(widget, setter)
                try:
                    fn(value)  # Set the widget.
                except Exception as e:
                    print(e) # handle type error

    def update_settings_from_widgets(self, map):
        for name, widget in map.items():
            cls = widget.__class__.__name__
            getter, setter, dtype = self.widget_mappers.get(cls, (None, None))
            print("save:", getter, setter)
            if getter:
                fn = getattr(widget, getter)
                value = fn()                
                print("-- value:", value, type(value), dtype)
                if value is not None:
                    self.settings.setValue(name, value) # Set the settings.

        # Notify watcher of changed settings.
        self.settings_changed.emit()
