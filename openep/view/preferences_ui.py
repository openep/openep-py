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

from unittest.util import sorted_list_difference
from PySide2 import QtWidgets, QtCore, QtGui
from PySide2.QtCore import Qt
from numpy import sign
import qdarkstyle

from openep.view.custom_widgets import CustomDockWidget


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
        self.tabs = QtWidgets.QTabWidget()
        self.tabs.setTabPosition(self.tabs.North)
        self.create_3d_viewer_settings()
        self.create_table_settings()
        self.create_annotation_settings()
        self.create_interpolation_settings()

        # TODO:
        # Add tab for the mapping points tables.
        # Default columns visible
        # Default sorting (column name, ascending or descinding)
        # Interpolate on deleting/restoring mapping points

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
        select_none  = QtWidgets.QRadioButton("Off")
        select_3d.setChecked(True)
        select.addButton(select_3d)
        select.addButton(select_surface)
        select.addButton(select_none)
        select.setExclusive(True)

        point_selection = QtWidgets.QGroupBox("Point selection")
        select_row = QtWidgets.QHBoxLayout()
        select_text_column = QtWidgets.QVBoxLayout()
        select_text_column.setAlignment(Qt.AlignTop)
        select_text = QtWidgets.QLabel("Select mapping points using the left mouse button:")
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

    def create_table_settings(self):
        """Settings for the mapping points and recycle bin tables."""

        # Mapping points table
        index = QtWidgets.QCheckBox("Index")
        tag = QtWidgets.QCheckBox("Tag")
        name = QtWidgets.QCheckBox("Name")
        voltage = QtWidgets.QCheckBox("Voltage")
        lat = QtWidgets.QCheckBox("LAT")
        boxes = [index, tag, name, voltage, lat]

        mapping_points = QtWidgets.QGroupBox("Mapping points")
        mapping_points_row = QtWidgets.QHBoxLayout()
        mapping_points_text_column = QtWidgets.QVBoxLayout()
        mapping_points_text_column.setAlignment(Qt.AlignTop)
        mapping_points_text = QtWidgets.QLabel("Columns visible at startup:")
        mapping_points_text_column.addWidget(mapping_points_text)
        mapping_points_row.addLayout(mapping_points_text_column, 0)
        mapping_points_checkboxes_layout = QtWidgets.QVBoxLayout()
        for box in boxes:
            box.setChecked(True)
            mapping_points_checkboxes_layout.addWidget(box)
        mapping_points_row.addLayout(mapping_points_checkboxes_layout, 0)
        mapping_points_row.addStretch(1)
        mapping_points.setLayout(mapping_points_row)

        # Recycle bin table
        index = QtWidgets.QCheckBox("Index")
        tag = QtWidgets.QCheckBox("Tag")
        name = QtWidgets.QCheckBox("Name")
        voltage = QtWidgets.QCheckBox("Voltage")
        lat = QtWidgets.QCheckBox("LAT")
        boxes = [index, tag, name, voltage, lat]

        recycle_bin = QtWidgets.QGroupBox("Recycle bin")
        recycle_bin_row = QtWidgets.QHBoxLayout()
        recycle_bin_text_column = QtWidgets.QVBoxLayout()
        recycle_bin_text_column.setAlignment(Qt.AlignTop)
        recycle_bin_text = QtWidgets.QLabel("Columns visible at startup:")
        recycle_bin_text_column.addWidget(recycle_bin_text)
        recycle_bin_row.addLayout(recycle_bin_text_column, 0)
        recycle_bin_checkboxes_layout = QtWidgets.QVBoxLayout()
        for box in boxes:
            box.setChecked(True)
            recycle_bin_checkboxes_layout.addWidget(box)
        recycle_bin_row.addLayout(recycle_bin_checkboxes_layout, 0)
        recycle_bin_row.addStretch(1)
        recycle_bin.setLayout(recycle_bin_row)

        # Sorting
        sorting = QtWidgets.QGroupBox("Default sorting")
        sorting_layout = QtWidgets.QHBoxLayout()

        sort_by_text = QtWidgets.QLabel("Sort by:")
        sort_by = QtWidgets.QComboBox()
        sort_by.addItems([
            "Index",
            "Tag",
            "Name",
            "Voltage",
            "LAT",
        ])
        sort_by.setMinimumWidth(120)
        # TODO: can use sort_by.setCurrentText('TEXT') to set the default selection
        
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

        # Set nested layouts
        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(mapping_points)
        layout.addWidget(recycle_bin)
        layout.addWidget(sorting)
        layout.addStretch()

        viewer = QtWidgets.QWidget()
        viewer.setLayout(layout)
        self.tabs.addTab(viewer, "Tables")

    def create_annotation_settings(self):
        """Settings for using the annotation viewer."""

        layout = QtWidgets.QVBoxLayout()

        # TODO: set validators for QLineEdit
        # See: https://zetcode.com/pyqt/qlineedit/
        # or https://stackoverflow.com/a/13423244/17623640

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

        # TODO: Add QStackLayout for the Signals and Annotaitons options
        #       Add QDoubleSpinbox for line thicknesses.
        #       Signals: Active and Non-active line thickness
        #       Annotations: WOI line thickness
        #                   Annotation Line thickness
        #                   Annotation size
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

        # Linewidth: Active: [] Others: []
        # Linewidth: [] Point size: 
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
        gain_options.addRow(QtWidgets.QLabel("Min."), QtWidgets.QLineEdit())
        gain_options.addRow(QtWidgets.QLabel("Max."), QtWidgets.QLineEdit())
        gain_options.addRow(QtWidgets.QLabel("Scroll speed."), QtWidgets.QLineEdit())
        gain_row.addLayout(gain_options, 0)
        gain_row.addStretch(1)
        gain.setLayout(gain_row)

        # Interpolation settings
        interpolate = QtWidgets.QCheckBox("Interpolate after re-annotating")
        interpolate.setCheckable(True)
        interpolate.setChecked(False)

        #layout.addWidget(figure)
        layout.addWidget(lines, 0)
        layout.addWidget(gain, 0)
        layout.addWidget(interpolate, 0)
        layout.addStretch(1)

        annotate = QtWidgets.QWidget()
        annotate.setLayout(layout)
        self.tabs.addTab(annotate, "Annotate")

    def create_interpolation_settings(self):
        """Settings for interpolation method and parameters."""

        method = QtWidgets.QGroupBox("Method")

        select_method = QtWidgets.QComboBox()
        select_method.setMinimumWidth(120)
        select_method.setEditable(False)
        select_method.addItem("RBF")

        select_method_layout = QtWidgets.QHBoxLayout()
        select_method_layout.addWidget(select_method)
        select_method_layout.addStretch()

        # TODO: Add options for RBF.
        # See: https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.RBFInterpolator.html
        method_options = QtWidgets.QGroupBox("Parameters")
        method_options_layout = QtWidgets.QStackedLayout()
        method_options.setLayout(method_options_layout)
        
        method_layout = QtWidgets.QVBoxLayout()
        method_layout.addLayout(select_method_layout)
        method.setLayout(method_layout)
        
        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(method, 0)  # do not expand vertical size
        layout.addWidget(method_options, 0)
        layout.addStretch(1)  # expand vertically to take up all remaining space

        interpolation = QtWidgets.QWidget()        
        interpolation.setLayout(layout)
        self.tabs.addTab(interpolation, "Interpolation")

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

    def enable_apply_or_discard(self):
        self.apply_or_discard.setEnabled(True)

class MainWindow(QtWidgets.QMainWindow):

    def __init__(self):

        super().__init__()
        self.preferences = PreferencesWidget("Preferences")
        self.addDockWidget(QtCore.Qt.RightDockWidgetArea, self.preferences)

        self.preferences.apply_or_discard.accepted.connect(self.accept)
        self.preferences.apply_or_discard.rejected.connect(self.reject)

    def accept(self):
        self.preferences.apply_or_discard.setEnabled(False)
        self.blockSignals(True)
        # TODO: Update the QSettings based on the values of the widgets.
        #       Also updateUpdate:
        #              * 3D viewers: colourmaps, colour and size of mapping points, interpolate if method or parameters changed
        #              * annotater: line thicknesses, min and max gain

        # Copy setting values from widgets to the settings manager
        # self.user_settings.update_settings_from_widgets(self.map)

        # Write the settings to file
        # self.user_settings.sync()

        self.blockSignals(False)

    def reject(self):
        self.preferences.apply_or_discard.setEnabled(False)
        self.blockSignals(True)
        # TODO: Update the values of the widgets based on QSettings.
        
        # Push previous settings back to the widgets
        # self.user_settings.update_widgets_from_settings(self.map)

        self.blockSignals(False)

if __name__ == '__main__':

    app = QtWidgets.QApplication([])

    app.setStyle("Fusion")

    # setup stylesheet
    app.setStyleSheet(qdarkstyle.load_stylesheet(qt_api='pyside2'))

    w = MainWindow()
    w.showMaximized()

    app.exec_()
