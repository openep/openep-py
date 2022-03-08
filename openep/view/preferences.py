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

from PySide2 import QtCore

__all__ = ["PreferencesManager"]


class PreferencesManager(QtCore.QSettings):
    """Settings manager for OpenEP Gui."""

    widget_mappers = {
        'QCheckBox': ('isChecked', 'setChecked', bool),
        'QLineEdit': ('text', 'setText', str),
        'QSpinBox': ('value', 'setValue', int),
        'QRadioButton': ('isChecked', 'setChecked', bool),
        'QComboBox': ('currentText', 'setCurrentText', str),
        'QSpinBox': ('value', 'setValue', int),
        'QDoubleSpinBox': ('value', 'setValue', float),
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
            print(self.widget_mappers.get(cls, (None, None)), cls)
            getter, setter, dtype = self.widget_mappers.get(cls, (None, None))
            value = self.settings.value(name, type=dtype)
            print("load:", getter, setter, value, type(value), dtype)
            if setter and value is not None:
                fn = getattr(widget, setter)
                try:
                    fn(value)  # Set the widget.
                except Exception as e:
                    print(e)  # handle type error

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
                    self.settings.setValue(name, value)  # Set the settings.

        # Notify watcher of changed settings.
        self.settings_changed.emit()

    def extract_preferences(self):
        """Create dictionary of preferences."""

        data = {}

        # 3D viewers settings
        # 3D viewers settings: Point selection
        data['3DViewers/PointSelection'] = self.settings.value('3DViewers/PointSelection')
        data['3DViewers/PointSelection/3D'] = self.settings.value('3DViewers/PointSelection/3D')
        data['3DViewers/PointSelection/Surface'] = self.settings.value('3DViewers/PointSelection/Surface')
        data['3DViewers/PointSelection/Off'] = self.settings.value('3DViewers/PointSelection/Off')

        point_selection_id = self.settings.value('3DViewers/PointSelection')
        print(point_selection_id)
        if point_selection_id == 0:
            data['3DViewers/PointSelection'] = "3DPoints"
        elif point_selection_id == 1:
            data['3DViewers/PointSelection'] = "ProjectedPoints"
        elif point_selection_id == 2:
            data['3DViewers/PointSelection'] = "Off"

        # 3D viewers settings: Secondary viewers
        link_secondary_viewers = self.settings.value('3DViewers/SecondaryViewers/Link')
        data['3DViewers/SecondaryViewers/Link'] = link_secondary_viewers

        # Table settings
        columns = ['Index', 'Tag', 'Name', 'Voltage', 'LAT']

        # Table settings: Mapping points
        data[f'Tables/MappingPoints/Hide'] = []
        for column in columns:
            show = self.settings.value(f'Tables/MappingPoints/Show/{column}')
            if not show:
                data[f'Tables/MappingPoints/Hide'].append(column)

        # Table settings: Recycle bin
        data[f'Tables/RecycleBin/Hide'] = []
        for column in columns:
            show = self.settings.value(f'Tables/RecycleBin/Show/{column}')
            if not show:
                data[f'Tables/RecycleBin/Hide'].append(column)

        # Table settings: Sorting
        sort_by = self.settings.value('Tables/SortBy')
        sort_order_text = self.settings.value('Tables/SortOrder')
        sort_order = QtCore.Qt.AscendingOrder if sort_order_text == "Ascending" else QtCore.Qt.DescendingOrder
        data['Tables/SortBy'] = sort_by
        data['Tables/SortOrder'] = sort_order

        # Table settings: Interpolate
        interpolate = self.settings.value('Tables/Interpolate')
        data['Tables/Interpolate'] = interpolate

        # Annotation settings
        data['Annotate/Lines/Signals/Linewidth/Active'] = self.settings.value('Annotate/Lines/Signals/Linewidth/Active')
        data['Annotate/Lines/Signals/Linewidth/Other'] = self.settings.value('Annotate/Lines/Signals/Linewidth/Other')
        data['Annotate/Lines/Annotations/Linewidth'] = self.settings.value('Annotate/Lines/Annotations/Linewidth')
        data['Annotate/Lines/Annotations/Markersize'] = self.settings.value('Annotate/Lines/Annotations/Markersize')

        data['Annotate/Gain/Min'] = float(self.settings.value('Annotate/Gain/Min'))
        data['Annotate/Gain/Max'] = float(self.settings.value('Annotate/Gain/Max'))
        data['Annotate/Gain/Prefactor'] = float(self.settings.value('Annotate/Gain/Prefactor'))

        data['Annotate/Interpolate'] = self.settings.value('Annotate/Interpolate')

        # Interpolation settings
        data['Interpolation/Method'] = self.settings.value('Interpolation/Method')
        data['Interpolation/RBFParameters/Neighbours'] = int(self.settings.value('Interpolation/RBFParameters/Neighbours'))
        if data['Interpolation/RBFParameters/Neighbours'] == 0:
            data['Interpolation/RBFParameters/Neighbours'] = None
        data['Interpolation/RBFParameters/Smoothing'] = float(self.settings.value('Interpolation/RBFParameters/Smoothing'))
        data['Interpolation/RBFParameters/Kernel'] = self.settings.value('Interpolation/RBFParameters/Kernel')
        data['Interpolation/RBFParameters/Epsilon'] = float(self.settings.value('Interpolation/RBFParameters/Epsilon'))
        data['Interpolation/RBFParameters/Degree'] = int(self.settings.value('Interpolation/RBFParameters/Degree'))

        return data
