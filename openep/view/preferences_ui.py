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
Create a dock widget for the setting preferences.
"""

from PySide2 import QtWidgets, QtGui, QtCore
from PySide2.QtCore import Qt
from scipy.interpolate import RBFInterpolator

from openep.view.custom_widgets import CustomDockWidget

__all__ = ['PreferencesWidget']

class PreferencesWidget(CustomDockWidget):

    def __init__(self, title: str):

        super().__init__(title)
        self.setMinimumSize(200, 100)

        # Make the preferences widget scrollable
        self.scroll = QtWidgets.QScrollArea()
        self.set_scroll_properties()

        # Create bold font for highlighting some text
        self.bold_font = QtGui.QFont()
        self.bold_font.setBold(True)

        # Add widgets for settings
        self.map = {}  # store all settings widgets in here so we can modify/retrieve values
        self.tabs = QtWidgets.QTabWidget()
        self.tabs.setTabPosition(self.tabs.North)
        self.create_3d_viewer_settings()
        self.create_table_settings()
        self.create_annotation_settings()
        self.create_interpolation_settings()

        # Add widget for accepting/rejecting changes to the settings
        self.create_accept_discard_buttons()

        # Set the overall layout
        self.set_nested_layout()

        # When any widget is modified, enable the Apply and Discard
        self.connect_all_signals()

    def set_scroll_properties(self):
        """Make the widget scrollable. Both vertically and horizontally."""

        self.scroll.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        self.scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        self.scroll.setWidgetResizable(True)

    def create_3d_viewer_settings(self):
        """Settings for the 3d viewers. Primary and secondary viewers have some distinct settings."""

        select = QtWidgets.QButtonGroup()
        select_3d = QtWidgets.QRadioButton("3D points")
        select_surface = QtWidgets.QRadioButton("Projected points")
        select_none = QtWidgets.QRadioButton("Off")
        select_3d.setChecked(True)
        select.addButton(select_3d)
        select.addButton(select_surface)
        select.addButton(select_none)
        select.setId(select_3d, 0)
        select.setId(select_surface, 1)
        select.setId(select_none, 2)
        select.setExclusive(True)

        point_selection = QtWidgets.QGroupBox("Point selection")
        select_row = QtWidgets.QHBoxLayout()
        select_text_column = QtWidgets.QVBoxLayout()
        select_text_column.setAlignment(Qt.AlignTop)
        select_text = QtWidgets.QLabel("Select mapping points using the left mouse button")
        select_text_column.addWidget(select_text)
        select_row.addLayout(select_text_column, 0)
        select_buttons_layout = QtWidgets.QVBoxLayout()
        select_buttons_layout.addWidget(select_3d)
        select_buttons_layout.addWidget(select_surface)
        select_buttons_layout.addWidget(select_none)
        select_row.addLayout(select_buttons_layout, 0)
        select_row.addStretch(1)
        point_selection.setLayout(select_row)

        secondary_viewers = QtWidgets.QGroupBox("Secondary viewers")
        secondary_viewers_layout = QtWidgets.QVBoxLayout()
        link = QtWidgets.QCheckBox("Link to primary viewer by default")
        link.setCheckable(True)
        link.setChecked(True)
        secondary_viewers_layout.addWidget(link)
        secondary_viewers.setLayout(secondary_viewers_layout)

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(point_selection)
        layout.addWidget(secondary_viewers)
        layout.addStretch()

        viewer = QtWidgets.QWidget()
        viewer.setLayout(layout)
        self.tabs.addTab(viewer, "3D viewers")

        # Add widgets to map
        self.map['3DViewers/PointSelection'] = select
        self.map['3DViewers/PointSelection/3D'] = select_3d
        self.map['3DViewers/PointSelection/Surface'] = select_surface
        self.map['3DViewers/PointSelection/Off'] = select_none
        self.map['3DViewers/SecondaryViewers/Link'] = link

    def create_table_settings(self):
        """Settings for the mapping points and recycle bin tables."""

        # Mapping points table
        mapping_points_index = QtWidgets.QCheckBox("Index")
        mapping_points_tag = QtWidgets.QCheckBox("Tag")
        mapping_points_name = QtWidgets.QCheckBox("Name")
        mapping_points_voltage = QtWidgets.QCheckBox("Voltage")
        mapping_points_lat = QtWidgets.QCheckBox("LAT")
        mapping_points_boxes = [
            mapping_points_index,
            mapping_points_tag,
            mapping_points_name,
            mapping_points_voltage,
            mapping_points_lat,
        ]

        mapping_points = QtWidgets.QGroupBox("Mapping points")
        mapping_points_row = QtWidgets.QHBoxLayout()
        mapping_points_text_column = QtWidgets.QVBoxLayout()
        mapping_points_text_column.setAlignment(Qt.AlignTop)
        mapping_points_text = QtWidgets.QLabel("Columns visible at startup")
        mapping_points_text_column.addWidget(mapping_points_text)
        mapping_points_row.addLayout(mapping_points_text_column, 0)
        mapping_points_checkboxes_layout = QtWidgets.QVBoxLayout()
        for box in mapping_points_boxes:
            box.setChecked(True)
            mapping_points_checkboxes_layout.addWidget(box)
        mapping_points_row.addLayout(mapping_points_checkboxes_layout, 0)
        mapping_points_row.addStretch(1)
        mapping_points.setLayout(mapping_points_row)

        # Recycle bin table
        recycle_bin_index = QtWidgets.QCheckBox("Index")
        recycle_bin_tag = QtWidgets.QCheckBox("Tag")
        recycle_bin_name = QtWidgets.QCheckBox("Name")
        recycle_bin_voltage = QtWidgets.QCheckBox("Voltage")
        recycle_bin_lat = QtWidgets.QCheckBox("LAT")
        recycle_bin_boxes = [
            recycle_bin_index,
            recycle_bin_tag,
            recycle_bin_name,
            recycle_bin_voltage,
            recycle_bin_lat,
        ]

        recycle_bin = QtWidgets.QGroupBox("Recycle bin")
        recycle_bin_row = QtWidgets.QHBoxLayout()
        recycle_bin_text_column = QtWidgets.QVBoxLayout()
        recycle_bin_text_column.setAlignment(Qt.AlignTop)
        recycle_bin_text = QtWidgets.QLabel("Columns visible at startup")
        recycle_bin_text_column.addWidget(recycle_bin_text)
        recycle_bin_row.addLayout(recycle_bin_text_column, 0)
        recycle_bin_checkboxes_layout = QtWidgets.QVBoxLayout()
        for box in recycle_bin_boxes:
            box.setChecked(True)
            recycle_bin_checkboxes_layout.addWidget(box)
        recycle_bin_row.addLayout(recycle_bin_checkboxes_layout, 0)
        recycle_bin_row.addStretch(1)
        recycle_bin.setLayout(recycle_bin_row)

        # Sorting
        sorting = QtWidgets.QGroupBox("Default sorting")
        sorting_layout = QtWidgets.QHBoxLayout()

        sort_by_text = QtWidgets.QLabel("Sort by")
        sort_by = QtWidgets.QComboBox()
        sort_by.addItems([
            "Index",
            "Tag",
            "Name",
            "Voltage",
            "LAT",
        ])
        sort_by.setMinimumWidth(120)

        sort_order_text = QtWidgets.QLabel("Order")
        sort_order = QtWidgets.QComboBox()
        sort_order.addItems(["Ascending", "Descending"])
        sort_order.setMinimumWidth(120)

        sorting_layout.addWidget(sort_by_text, 0)
        sorting_layout.addWidget(sort_by, 0)
        sorting_layout.addWidget(sort_order_text, 0)
        sorting_layout.addWidget(sort_order, 0)
        sorting_layout.addStretch(1)
        sorting.setLayout(sorting_layout)

        # Interpolation settings
        interpolate = QtWidgets.QCheckBox("Interpolate after deleting or restoring points.")
        interpolate.setCheckable(True)
        interpolate.setChecked(True)

        # Set nested layouts
        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(mapping_points)
        layout.addWidget(recycle_bin)
        layout.addWidget(sorting)
        layout.addWidget(interpolate)
        layout.addStretch()

        viewer = QtWidgets.QWidget()
        viewer.setLayout(layout)
        self.tabs.addTab(viewer, "Tables")

        # Add widgets to map
        self.map['Tables/MappingPoints/Show/Index'] = mapping_points_index
        self.map['Tables/MappingPoints/Show/Tag'] = mapping_points_tag
        self.map['Tables/MappingPoints/Show/Name'] = mapping_points_name
        self.map['Tables/MappingPoints/Show/Voltage'] = mapping_points_voltage
        self.map['Tables/MappingPoints/Show/LAT'] = mapping_points_lat

        self.map['Tables/RecycleBin/Show/Index'] = recycle_bin_index
        self.map['Tables/RecycleBin/Show/Tag'] = recycle_bin_tag
        self.map['Tables/RecycleBin/Show/Name'] = recycle_bin_name
        self.map['Tables/RecycleBin/Show/Voltage'] = recycle_bin_voltage
        self.map['Tables/RecycleBin/Show/LAT'] = recycle_bin_lat

        self.map['Tables/SortBy'] = sort_by
        self.map['Tables/SortOrder'] = sort_order

        self.map['Tables/Interpolate'] = interpolate

    def create_annotation_settings(self):
        """Settings for using the annotation viewer."""

        layout = QtWidgets.QVBoxLayout()

        # TODO: Add legend section
        #       Fontsize
        #       Linewidth
        #       Draggable (move the legend around the canvas), False

        # TODO: Add these options, provide recommended settings based on screen resolution
        figure = QtWidgets.QGroupBox("Figure")
        figure_row = QtWidgets.QHBoxLayout()
        figure_options = QtWidgets.QFormLayout()
        figure_options.addRow(QtWidgets.QLabel("Height"), QtWidgets.QLineEdit())
        figure_options.addRow(QtWidgets.QLabel("Width"), QtWidgets.QLineEdit())
        figure_options.addRow(QtWidgets.QLabel("DPI"), QtWidgets.QLineEdit())
        figure_row.addLayout(figure_options, 0)
        figure_row.addStretch(1)
        figure.setLayout(figure_row)

        # Linestyle options
        lines = QtWidgets.QGroupBox("Lines")
        lines_layout = QtWidgets.QVBoxLayout()

        lines_row = QtWidgets.QHBoxLayout()
        line_types = QtWidgets.QComboBox()
        line_types.addItems(["Signals", "Annotations"])
        line_types.setMinimumWidth(120)
        lines_row.addWidget(line_types)
        lines_row.addStretch(1)

        line_options = QtWidgets.QStackedLayout()
        line_types.activated.connect(line_options.setCurrentIndex)

        # Linestyle options: Signals
        signal_options = QtWidgets.QWidget()
        signal_options_layout = QtWidgets.QVBoxLayout()
        signal_linewidth_heading = QtWidgets.QLabel("Linewidth")
        signal_linewidth_heading.setFont(self.bold_font)
        signal_options_layout.addWidget(signal_linewidth_heading)

        signal_linewidth_layout = QtWidgets.QHBoxLayout()
        active_linewidth_text = QtWidgets.QLabel("Active: ")
        active_linewidth = QtWidgets.QDoubleSpinBox()
        active_linewidth.setFixedWidth(60)
        active_linewidth.setRange(0.5, 5.0)
        active_linewidth.setSingleStep(0.1)
        active_linewidth.setDecimals(1)
        active_linewidth.setWrapping(False)
        active_linewidth.setValue(1.2)
        active_linewidth.textFromValue(1.2)
        other_linewidth_text = QtWidgets.QLabel("Others: ")
        other_linewidth = QtWidgets.QDoubleSpinBox()
        other_linewidth.setFixedWidth(60)
        other_linewidth.setRange(0.5, 5.0)
        other_linewidth.setSingleStep(0.1)
        other_linewidth.setDecimals(1)
        other_linewidth.setWrapping(False)
        other_linewidth.setValue(0.8)
        other_linewidth.textFromValue(0.8)

        signal_linewidth_layout.addWidget(active_linewidth_text, 0)
        signal_linewidth_layout.addWidget(active_linewidth, 0)
        signal_linewidth_layout.addWidget(other_linewidth_text, 0)
        signal_linewidth_layout.addWidget(other_linewidth, 0)
        signal_linewidth_layout.addStretch(1)
        signal_options_layout.addLayout(signal_linewidth_layout, 0)
        signal_options_layout.addStretch(1)
        signal_options.setLayout(signal_options_layout)
        line_options.addWidget(signal_options)

        # Linestyle options: Annotations
        annotation_options = QtWidgets.QWidget()
        annotation_options_layout = QtWidgets.QVBoxLayout()
        annotation_linewidth_heading = QtWidgets.QLabel("Linewidth and size")
        annotation_linewidth_heading.setFont(self.bold_font)
        annotation_options_layout.addWidget(annotation_linewidth_heading)

        annotation_layout = QtWidgets.QHBoxLayout()
        annotation_linewidth_text = QtWidgets.QLabel("Linewidth: ")
        annotation_linewidth = QtWidgets.QDoubleSpinBox()
        annotation_linewidth.setFixedWidth(60)
        annotation_linewidth.setRange(0.5, 5.0)
        annotation_linewidth.setSingleStep(0.1)
        annotation_linewidth.setDecimals(1)
        annotation_linewidth.setWrapping(False)
        annotation_linewidth.setValue(1.0)
        annotation_linewidth.textFromValue(1.0)
        annotation_size_text = QtWidgets.QLabel("Marker size: ")
        annotation_size = QtWidgets.QDoubleSpinBox()
        annotation_size.setFixedWidth(60)
        annotation_size.setRange(1, 10)
        annotation_size.setSingleStep(0.2)
        annotation_size.setDecimals(1)
        annotation_size.setWrapping(False)
        annotation_size.setValue(4)
        annotation_size.textFromValue(4)

        annotation_layout.addWidget(annotation_linewidth_text, 0)
        annotation_layout.addWidget(annotation_linewidth, 0)
        annotation_layout.addWidget(annotation_size_text, 0)
        annotation_layout.addWidget(annotation_size, 0)
        annotation_layout.addStretch(1)
        annotation_options_layout.addLayout(annotation_layout, 0)
        annotation_options_layout.addStretch(1)
        annotation_options.setLayout(annotation_options_layout)
        line_options.addWidget(annotation_options)

        # Linestyle options: Set nested layout
        lines_layout.addLayout(lines_row, 0)
        lines_layout.addLayout(line_options, 0)
        lines_layout.addStretch(1)
        lines.setLayout(lines_layout)

        # Signal gain
        gain = QtWidgets.QGroupBox("Gain")
        gain_row = QtWidgets.QHBoxLayout()
        gain_options = QtWidgets.QFormLayout()

        min_gain = QtWidgets.QLineEdit("-5.0")
        min_gain_validator = QtGui.QDoubleValidator()
        min_gain_validator.setRange(-1e6, 1e6)
        min_gain_validator.setDecimals(1)
        min_gain.setValidator(min_gain_validator)
        gain_options.addRow(QtWidgets.QLabel("Min."), min_gain)

        max_gain = QtWidgets.QLineEdit("2.0")
        max_gain_validator = QtGui.QDoubleValidator()
        max_gain_validator.setRange(-1e6, 1e6)
        max_gain_validator.setDecimals(1)
        max_gain.setValidator(max_gain_validator)
        gain_options.addRow(QtWidgets.QLabel("Max."), max_gain)

        scroll_speed = QtWidgets.QLineEdit("0.001")
        scroll_speed_validator = QtGui.QDoubleValidator()
        scroll_speed_validator.setRange(0.01, 1)
        scroll_speed_validator.setDecimals(3)
        scroll_speed.setValidator(scroll_speed_validator)
        gain_options.addRow(QtWidgets.QLabel("Scroll speed."), scroll_speed)

        gain_row.addLayout(gain_options, 0)
        gain_row.addStretch(1)
        gain.setLayout(gain_row)

        # Interpolation settings
        interpolate = QtWidgets.QCheckBox("Interpolate after re-annotating")
        interpolate.setCheckable(True)
        interpolate.setChecked(False)

        # layout.addWidget(figure)
        layout.addWidget(lines, 0)
        layout.addWidget(gain, 0)
        layout.addWidget(interpolate, 0)
        layout.addStretch(1)

        annotate = QtWidgets.QWidget()
        annotate.setLayout(layout)
        self.tabs.addTab(annotate, "Annotate")

        # Add widgets to map
        self.map['Annotate/Lines/Signals/Linewidth/Active'] = active_linewidth
        self.map['Annotate/Lines/Signals/Linewidth/Other'] = other_linewidth
        self.map['Annotate/Lines/Annotations/Linewidth'] = annotation_linewidth
        self.map['Annotate/Lines/Annotations/Markersize'] = annotation_size

        self.map['Annotate/Gain/Min'] = min_gain
        self.map['Annotate/Gain/Max'] = max_gain
        self.map['Annotate/Gain/Prefactor'] = scroll_speed

        self.map['Annotate/Interpolate'] = interpolate

    def create_interpolation_settings(self):
        """Settings for interpolation method and parameters."""

        # Selecting the interpolation method
        method = QtWidgets.QGroupBox("Method")
        select_method = QtWidgets.QComboBox()
        select_method.setMinimumWidth(120)
        select_method.setEditable(False)
        select_method.addItem("RBF")

        select_method_layout = QtWidgets.QHBoxLayout()
        select_method_layout.addWidget(select_method)
        select_method_layout.addStretch()

        method_layout = QtWidgets.QVBoxLayout()
        method_layout.addLayout(select_method_layout)
        method.setLayout(method_layout)

        # Setting the interpolation parameters
        method_options = QtWidgets.QGroupBox("Parameters")
        method_options_layout = QtWidgets.QStackedLayout()
        select_method.activated.connect(method_options_layout.setCurrentIndex)

        # Setting the interpolation parameters: RBF
        rbf = QtWidgets.QWidget()
        rbf_layout = QtWidgets.QVBoxLayout()
        rbf_parameters_layout = QtWidgets.QFormLayout()

        neighbours = QtWidgets.QLineEdit()
        neighbours.setText("0")
        neighbours_validator = QtGui.QIntValidator()
        neighbours_validator.setRange(0, 1e7)
        neighbours.setValidator(neighbours_validator)
        rbf_parameters_layout.addRow("Neighbours<sup>*</sup>", neighbours)

        smoothing = QtWidgets.QDoubleSpinBox()
        smoothing.setRange(0, 20)
        smoothing.setSingleStep(0.1)
        smoothing.setDecimals(1)
        smoothing.setWrapping(False)
        smoothing.setValue(4)
        rbf_parameters_layout.addRow("Smoothing", smoothing)

        kernel = QtWidgets.QComboBox()
        kernel.addItems([
            "linear",
            "thin_plate_spline",
            "cubic",
            "quintic",
            "multiquadric",
            "inverse_multiquadric",
            "inverse_quadratic",
            "gaussian",
        ])
        kernel.setCurrentText("multiquadric")
        rbf_parameters_layout.addRow("Kernel", kernel)

        epsilon = QtWidgets.QDoubleSpinBox()
        epsilon.setRange(0.1, 20)
        epsilon.setSingleStep(0.1)
        epsilon.setDecimals(1)
        epsilon.setWrapping(False)
        epsilon.setValue(1)
        rbf_parameters_layout.addRow("Epsilon", epsilon)

        # TODO: Add custom validator. Need to ensure degree is not too small
        #       See: https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.RBFInterpolator.html
        #       And: https://doc.qt.io/qt-5/qvalidator.html
        degree = QtWidgets.QSpinBox()
        degree.setRange(-1, 10)
        degree.setSingleStep(1)
        degree.setWrapping(False)
        degree.setValue(1)
        rbf_parameters_layout.addRow("Degree", degree)

        zero_is_none = QtWidgets.QLabel("<sup>*</sup>A value of 0 will include all neighbours.")
        scipy_docs = QtWidgets.QLabel()
        scipy_docs.setOpenExternalLinks(True)
        scipy_url = "https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.RBFInterpolator.html"
        scipy_text = "Scipy documenation"
        scipy_docs.setText(f"See the <a href={scipy_url}>{scipy_text}</a> for a description of the parameters.")

        rbf_parameters_row = QtWidgets.QHBoxLayout()
        rbf_parameters_row.addLayout(rbf_parameters_layout, 0)
        rbf_parameters_row.addStretch(1)

        rbf_layout.addLayout(rbf_parameters_row)
        rbf_layout.addWidget(zero_is_none)
        rbf_layout.addWidget(scipy_docs)
        rbf.setLayout(rbf_layout)

        # Add all parameter layouts
        method_options_layout.addWidget(rbf)

        # For some reason, the QStackedLayout overlaps with the QGroupBox title
        # Need to put the stacked layout inside a VBox so there's no overlap
        method_options_overall_layout = QtWidgets.QVBoxLayout()
        method_options_overall_layout.addLayout(method_options_layout)
        method_options.setLayout(method_options_overall_layout)

        # Setting the nested layout
        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(method, 0)  # do not expand vertical size
        layout.addWidget(method_options, 0)
        layout.addStretch(1)  # expand vertically to take up all remaining space

        interpolation = QtWidgets.QWidget()
        interpolation.setLayout(layout)
        self.tabs.addTab(interpolation, "Interpolation")

        # Add widgets to map
        self.map["Interpolation/Method"] = select_method
        self.map["Interpolation/RBFParameters/Neighbours"] = neighbours
        self.map["Interpolation/RBFParameters/Smoothing"] = smoothing
        self.map["Interpolation/RBFParameters/Kernel"] = kernel
        self.map["Interpolation/RBFParameters/Epsilon"] = epsilon
        self.map["Interpolation/RBFParameters/Degree"] = degree

    def create_accept_discard_buttons(self):
        """Buttons for applying or discarding the changes made to the settings."""

        _buttons = QtWidgets.QDialogButtonBox
        self.apply_or_discard = QtWidgets.QDialogButtonBox()
        apply = QtWidgets.QPushButton("Apply")
        discard = QtWidgets.QPushButton("Discard")
        self.apply_or_discard.addButton(apply, _buttons.AcceptRole)
        self.apply_or_discard.addButton(discard, _buttons.RejectRole)
        self.apply_or_discard.setEnabled(False)

    def set_nested_layout(self):
        """Create the layout for the Preferences dock."""

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.tabs)
        layout.addWidget(self.apply_or_discard)

        containing_widget = QtWidgets.QWidget()
        containing_widget.setLayout(layout)
        self.scroll.setWidget(containing_widget)
        self.setWidget(self.scroll)

    def connect_all_signals(self):
        """When the value of any widget is modified, enable the Apply and Discard buttons."""

        combo_boxes = self.tabs.findChildren(QtWidgets.QComboBox)
        for combo_box in combo_boxes:
            combo_box.currentIndexChanged.connect(self.enable_apply_or_discard)

        radio_buttons = self.tabs.findChildren(QtWidgets.QRadioButton)
        for radio_button in radio_buttons:
            radio_button.toggled.connect(self.enable_apply_or_discard)

        checkboxes = self.tabs.findChildren(QtWidgets.QCheckBox)
        for checkbox in checkboxes:
            checkbox.toggled.connect(self.enable_apply_or_discard)

        line_edits = self.tabs.findChildren(QtWidgets.QLineEdit)
        for line_edit in line_edits:
            line_edit.editingFinished.connect(self.enable_apply_or_discard)

        spin_boxes = self.tabs.findChildren(QtWidgets.QSpinBox)
        for spin_box in spin_boxes:
            spin_box.valueChanged.connect(self.enable_apply_or_discard)

        double_spin_boxes = self.tabs.findChildren(QtWidgets.QDoubleSpinBox)
        for double_spin_box in double_spin_boxes:
            double_spin_box.valueChanged.connect(self.enable_apply_or_discard)

    def enable_apply_or_discard(self):
        self.apply_or_discard.setEnabled(True)

    def extract_preferences(self):
        """Create dictionary of preferences."""

        data = {}

        # 3D viewers settings
        # 3D viewers settings: Point selection
        data['3DViewers/PointSelection/3D'] = self.map['3DViewers/PointSelection/3D'].isChecked()
        data['3DViewers/PointSelection/Surface'] = self.map['3DViewers/PointSelection/Surface'].isChecked()
        data['3DViewers/PointSelection/Off'] = self.map['3DViewers/PointSelection/Off'].isChecked()

        point_selection_id = self.map['3DViewers/PointSelection'].checkedId()
        if point_selection_id == 0:
            data['3DViewers/PointSelection'] = "3DPoints"
        elif point_selection_id == 1:
            data['3DViewers/PointSelection'] = "ProjectedPoints"
        elif point_selection_id == 2:
            data['3DViewers/PointSelection'] = "Off"

        # 3D viewers settings: Secondary viewers
        link_secondary_viewers = self.map['3DViewers/SecondaryViewers/Link'].isChecked()
        data['3DViewers/SecondaryViewers/Link'] = link_secondary_viewers

        # Table settings
        columns = ['Index', 'Tag', 'Name', 'Voltage', 'LAT']

        # Table settings: Mapping points
        data[f'Tables/MappingPoints/Hide'] = []
        for column in columns:
            show = self.map[f'Tables/MappingPoints/Show/{column}'].isChecked()
            if not show:
                data[f'Tables/MappingPoints/Hide'].append(column)

        # Table settings: Recycle bin
        data[f'Tables/RecycleBin/Hide'] = []
        for column in columns:
            show = self.map[f'Tables/RecycleBin/Show/{column}'].isChecked()
            if not show:
                data[f'Tables/RecycleBin/Hide'].append(column)

        # Table settings: Sorting
        sort_by = self.map['Tables/SortBy'].currentText()
        sort_order_text = self.map['Tables/SortOrder'].currentText()
        sort_order = QtCore.Qt.AscendingOrder if sort_order_text == "Ascending" else QtCore.Qt.DescendingOrder
        data['Tables/SortBy'] = sort_by
        data['Tables/SortOrder'] = sort_order

        # Table settings: Interpolate
        interpolate =self.map['Tables/Interpolate'].isChecked()
        data['Tables/Interpolate'] = interpolate

        # Annotation settings
        data['Annotate/Lines/Signals/Linewidth/Active'] = self.map['Annotate/Lines/Signals/Linewidth/Active'].value()
        data['Annotate/Lines/Signals/Linewidth/Other'] = self.map['Annotate/Lines/Signals/Linewidth/Other'].value()
        data['Annotate/Lines/Annotations/Linewidth'] = self.map['Annotate/Lines/Annotations/Linewidth'].value()
        data['Annotate/Lines/Annotations/Markersize'] = self.map['Annotate/Lines/Annotations/Markersize'].value()

        data['Annotate/Gain/Min'] = float(self.map['Annotate/Gain/Min'].text())
        data['Annotate/Gain/Max'] = float(self.map['Annotate/Gain/Max'].text())
        data['Annotate/Gain/Prefactor'] = float(self.map['Annotate/Gain/Prefactor'].text())

        self.map['Annotate/Interpolate'].isChecked()

        # Interpolation settings
        method_text = self.map['Interpolation/Method'].currentText()

        if method_text == "RBF":

            method = RBFInterpolator
            neighbours = int(self.map['Interpolation/RBFParameters/Neighbours'].text())
            neighbours = None if neighbours == 0 else neighbours
            parameters = {
                "neighbors": neighbours,
                "smoothing": self.map['Interpolation/RBFParameters/Smoothing'].value(),
                "kernel": self.map['Interpolation/RBFParameters/Kernel'].currentText(),
                "epsilon": self.map['Interpolation/RBFParameters/Epsilon'].value(),
                "degree": self.map['Interpolation/RBFParameters/Degree'].value(),
            }

        data['Interpolation/Method'] = method
        data['Interpolation/Parameters'] = parameters

        return data

