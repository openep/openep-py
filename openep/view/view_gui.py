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
A GUI for OpenEP-Py.
"""
import sys
import pathlib
import re

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import Qt
import numpy as np

import openep
import openep.view.custom_widgets
import openep.view.canvases
import openep.view.images
import openep.view.system_manager_ui
import openep.view.system_manager


class OpenEPGUI(QtWidgets.QMainWindow):

    def __init__(self):

        super().__init__()

        self._init_systems()
        self._init_ui()
        self._create_system_manager_dock()
        self._create_system_manager_menubar()
        self._create_system_manager_table()
        self._create_egm_canvas_dock()
        self._create_analysis_canvas_dock()
        self._add_dock_widgets()
        self._disable_dock_widgets()

    def _init_systems(self):
        """Containers and variables for keeping track of the systems loaded into the GUI."""

        # We need to keep track of all cases currently loaded (both OpenEP and openCARP)
        # Each system will keep track of its own docks/plotters/meshes
        self.systems = {}
        self._system_counter = 0  # we need a unique ID for each new system (IDs of deleted systems should not be reused)

        # We will use this to determine which egms to plot
        # All windows belonging to the main case have their menubars coloured blue
        self._active_system = None

        self.egm_slider = None
        self.egm_selection_filter = re.compile('\d+')  # noqa: W605

    def _init_ui(self):
        """
        Initialise the GUI.

        Lots of the layout is handled automatically by QtWidgets.QMainWindow.
        """

        self.setMinimumSize(800, 600)
        self.setWindowTitle("OpenEP: The open-source solution for electrophysiology data analysis")

    def _create_system_manager_dock(self):
        """Create a dockable widget that for managing systems loaded into the GUI."""

        self.system_dock = openep.view.custom_widgets.CustomDockWidget("Systems")
        self.system_main = QtWidgets.QMainWindow()
        self.system_main.setStyleSheet('QMainWindow {border: 0px; background-color: #d8dcd6;}')

        # The dock is set to have bold font (so the title stands out)
        # But all other widgets should have normal weighted font
        main_font = QtGui.QFont()
        main_font.setBold(False)
        self.system_main.setFont(main_font)

    def _create_system_manager_menubar(self):
        """
        Add a menu bar to the 'System' window of the GUI.

        The file menu has options for loading either an OpenEP dataset or an openCARP one, as well as
        loading auxillary data into an existing openCARP system.

        The view menu has options for creating a new 3D-viewer for an existing system.
        """

        # For loading new systems, deleting systems, adding a view of existing systems
        # Adding data to systems, exporting system to a different format
        self.system_main.menubar = self.system_main.menuBar()
        self.system_main.menubar.setNativeMenuBar(False)
        self.system_main.menubar.setStyleSheet(
            "QMenuBar{"
            "background-color: #95a3a6;"  # xkcd:cool grey
            "color: black;"
            "border: None;"
            "}"
        )

        file_menu = self.system_main.menubar.addMenu("File")

        load_menu = QtWidgets.QMenu("Load", self.system_main)
        load_openep_action = QtWidgets.QAction("OpenEP", self.system_main)
        load_openep_action.triggered.connect(self._load_openep)
        load_menu.addAction(load_openep_action)
        load_openCARP_action = QtWidgets.QAction("openCARP", self.system_main)
        load_openCARP_action.triggered.connect(self._load_openCARP)
        load_menu.addAction(load_openCARP_action)

        # The below will be used for adding data to openCARP datasets
        self.system_main.add_data_menu = QtWidgets.QMenu("Add data to", self.system_main)

        # The below will be used exporting OpenEP mesh data in openCARP format
        self.system_main.export_data_menu = QtWidgets.QMenu("Export", self.system_main)

        file_menu.addMenu(load_menu)
        file_menu.addMenu(self.system_main.add_data_menu)
        file_menu.addSeparator()
        file_menu.addMenu(self.system_main.export_data_menu)

        # The below will be used for creating a new 3D viewer for a given system
        view_menu = self.system_main.menubar.addMenu("View")
        self.system_main.add_view_menu = QtWidgets.QMenu("Add view for", self.system_main)
        view_menu.addMenu(self.system_main.add_view_menu)
        view_menu.addSeparator()

    def _create_system_manager_table(self):
        """
        Create a widget that will display information about each system loaded, and
        allow selection of an active system.
        """

        system_main_widget = QtWidgets.QWidget()
        self.system_main_layout = QtWidgets.QGridLayout(system_main_widget)

        # add a row containing heading labels to the grid
        openep.view.system_manager_ui.add_heading_bar(layout=self.system_main_layout)

        # we need a method for selecting the active system
        self._active_system_button_group = QtWidgets.QButtonGroup()

        # After adding all the rows we need to make sure they are at the top of the layout
        self.system_main.vertical_stretch = openep.view.system_manager_ui.create_vertical_stretch()
        self.system_main_layout.addItem(self.system_main.vertical_stretch)

        # Setting nested layouts
        self.system_main.setCentralWidget(system_main_widget)
        self.system_dock.setWidget(self.system_main)

    def update_system_manager_table(self, system):
        """Update the System manager table when a new system is loaded.

        We need to remove the vertical spacer, add a new row in the system table,
        then put the vertical spacer back.

        We also need to add the system to the 'File > Add data to' submenu and the
        'View > Add view for' submenu.
        """

        self.system_main_layout.removeItem(self.system_main.vertical_stretch)

        # Not sure why * 10 works. But if not multiplied by at least 10, there is vertical overlap between rows
        row_number = len(self.systems) * 10

        system.name_widget, *_, active_widget = openep.view.system_manager_ui.create_widgets_for_system_row(
            layout=self.system_main_layout,
            name=str(system.name),
            directory=str(system.folder),
            data_type=system.type,
            is_active=True if self._active_system.name == system.name else False,
            row_number=row_number,
            )

        system.name_widget.returnPressed.connect(lambda: self.update_system_name(system))
        active_widget.clicked.connect(lambda: self.change_active_system(system))
        self._active_system_button_group.addButton(active_widget)

        # Make sure the vertical spacing is correct again
        self.system_main.vertical_stretch = openep.view.system_manager_ui.create_vertical_stretch()
        self.system_main_layout.addItem(self.system_main.vertical_stretch)

    def update_system_name(self, system):
        """
        Change the name of a system.

        We also need to change the window title of its dock widgets.
        """

        new_name = system.name_widget.text()

        if new_name not in self.systems:
            self.systems[new_name] = self.systems.pop(system.name)
            self.systems[new_name].name = new_name
            for dock, mesh  in zip(system.docks, system.meshes):
                dock.setWindowTitle(f"{system.name}: {mesh.active_scalars_info.name}")
        else:

            system.name_widget.setText(system.name)

            error = QtWidgets.QMessageBox()
            error.setIcon(QtWidgets.QMessageBox.Critical)
            error.setText("Name Error")
            error.setInformativeText(
                "System names must be unique."
            )
            error.setWindowTitle("Error")
            error.exec_()

    def _create_egm_canvas_dock(self):
        """
        Create a dockable widget for plotting interactive electrograms with matplotlib.

        The user can select which points the electrograms will be plotted for, as well
        as the type(s) of electrograms to plot: reference, bipolar, unipolar A,
        unipolar B.

        The user can also change the window of interest with a RangeSlider,
        and setting the window of interest with a push button will re-interpolate
        the electrogram data onto the surface (for OpenEP datasets) or re-calculate
        voltages directly from electrograms (for openCARP datasets).
        """

        self.egm_dock = openep.view.custom_widgets.CustomDockWidget("Electrograms")
        self.egm_canvas, self.egm_figure, self.egm_axis = openep.view.canvases.create_canvas()
        egm_layout = QtWidgets.QVBoxLayout(self.egm_canvas)

        # Add widgets for manually selecting electrograms from point indices
        egm_selection_layout, self.egm_selection = openep.view.canvases.create_egm_selection_layout()
        self.egm_selection.returnPressed.connect(self.update_electrograms_from_text)
        egm_layout.addLayout(egm_selection_layout)

        # Add radio buttons to select bipolar, unipolar (A/B), and reference electrograms
        egm_type_layout, *checkboxes = openep.view.canvases.create_egm_type_layout()

        for box in checkboxes:
            box.stateChanged.connect(self.plot_electrograms)

        self.egm_reference_checkbox, self.egm_bipolar_checkbox, self.egm_unipolar_A_checkbox, self.egm_unipolar_B_checkbox = \
            checkboxes

        egm_layout.addLayout(egm_type_layout)
        egm_layout.addStretch()

        # Button for setting the window of interest
        self.set_woi_button = openep.view.canvases.add_woi_set_button(
            figure=self.egm_figure,
            axis=[0.782, 0.02, 0.09, 0.025],
        )
        self.set_woi_button.on_clicked(self.update_scalar_fields)

        # Create toolbar for saving, zooming etc.
        toolbar = openep.view.canvases.add_toolbar(
            canvas=self.egm_canvas,
            parent=self.egm_dock,
        )

        # Create a placeholder widget to hold our toolbar and canvas.
        canvas_widget = openep.view.canvases.create_canvas_widget(
            canvas=self.egm_canvas,
            toolbar=toolbar,
        )

        # Create a placeholder widget to hold our canvas, toolbar, and selection widgets
        # We're using a QMainWindow so we can easily add a menubar
        egm_canvas_main = QtWidgets.QMainWindow()
        egm_canvas_main.setCentralWidget(canvas_widget)

        # The dock is set to have bold font (so the title stands out)
        # But all other widgets should have normal weighted font
        main_font = QtGui.QFont()
        main_font.setBold(False)
        egm_canvas_main.setFont(main_font)

        self.egm_dock.setWidget(egm_canvas_main)

    def _create_analysis_canvas_dock(self):
        """
        Create a dockable widget for other matplotlib plots.

        For example, plotting a histogram of the surface area occupied by
        a range of bipolar voltages.
        """

        self.analysis_dock = openep.view.custom_widgets.CustomDockWidget("Analysis")
        self.analysis_canvas, self.analysis_figure, self.analysis_axis = openep.view.canvases.create_canvas()
        analysis_layout = QtWidgets.QVBoxLayout(self.analysis_canvas)

        analysis_layout.addStretch()

        # Create toolbar for saving, zooming etc.
        toolbar = openep.view.canvases.add_toolbar(
            canvas=self.analysis_canvas,
            parent=self.analysis_dock,
            keep_actions=['Save'],
        )

        # Create a placeholder widget to hold our toolbar and canvas.
        canvas_widget = openep.view.canvases.create_canvas_widget(
            canvas=self.analysis_canvas,
            toolbar=toolbar,
        )

        # Create a placeholder widget to hold our canvas, toolbar, and selection widgets
        analysis_canvas_main = QtWidgets.QMainWindow()
        analysis_canvas_main.setCentralWidget(canvas_widget)

        # The dock is set to have bold font (so the title stands out)
        # But all other widgets should have normal weighted font
        main_font = QtGui.QFont()
        main_font.setBold(False)
        analysis_canvas_main.setFont(main_font)

        self.analysis_dock.setWidget(analysis_canvas_main)

    def _add_dock_widgets(self):
        """
        Add dockable widgets to the main window.

        The two BackgroundPlotters are tabified, as are the two MPL canvases.
        """

        self.addDockWidget(Qt.LeftDockWidgetArea, self.egm_dock)

        self.addDockWidget(Qt.RightDockWidgetArea, self.analysis_dock)
        self.addDockWidget(Qt.RightDockWidgetArea, self.system_dock)
        self.tabifyDockWidget(self.analysis_dock, self.system_dock)

        for dock in [self.egm_dock, self.analysis_dock, self.system_dock]:
            dock.setAllowedAreas(Qt.AllDockWidgetAreas)

        self.setDockOptions(self.GroupedDragging | self.AllowTabbedDocks | self.AllowNestedDocks)
        self.setTabPosition(Qt.AllDockWidgetAreas, QtWidgets.QTabWidget.North)

    def _disable_dock_widgets(self):
        """
        Ignore all key presses in the dock widgets.

        This is required if there is no file loaded, otherwise
        the GUI will crash.
        """

        self.egm_dock.setEnabled(False)
        self.analysis_dock.setEnabled(False)

    def _load_openep(self):
        """
        Load an OpenEP case.

        Currently, only MATLAB files are supported.
        """

        dialogue = QtWidgets.QFileDialog()
        dialogue.setWindowTitle('Load an OpenEP file')
        dialogue.setDirectory(QtCore.QDir.currentPath())
        dialogue.setFileMode(QtWidgets.QFileDialog.ExistingFile)
        dialogue.setNameFilter("MATLAB file (*.mat)")

        if dialogue.exec_():
            filename = dialogue.selectedFiles()[0]
            self._initialise_case(filename)

    def _initialise_case(self, filename):
        """Initialise data from an OpenEP case object.

        Create separate meshes for the two BackgroundPlotters.
        Interpolate data ready to be plotted.
        """

        case = openep.load_openep(filename)

        new_system = openep.view.system_manager.System(
            name=str(self._system_counter),
            folder=pathlib.Path(filename).parent.resolve(),
            type="OpenEP",
            data=case,
        )

        # we need a unique key for each system
        while self._system_counter in self.systems:
            self._system_counter += 1

        self.systems[new_system.name] = new_system
        self._system_counter += 1
        self._active_system = new_system if self._active_system is None else self._active_system

        self.update_system_manager_table(new_system)

        # We need to dynamically add options for exporting data/creating 3D viewers to the main menubar
        add_export_system_menu = QtWidgets.QMenu(filename, self)
        add_export_openCARP_action = QtWidgets.QAction("as openCARP", self)
        add_export_openCARP_action.triggered.connect(lambda: self._export_data_to_openCARP(new_system))
        add_export_system_menu.addAction(add_export_openCARP_action)
        self.system_main.export_data_menu.addMenu(add_export_system_menu)

        add_view_for_system_action = QtWidgets.QAction(filename, self)
        add_view_for_system_action.triggered.connect(lambda: self.add_view(new_system))
        self.system_main.add_view_menu.addAction(add_view_for_system_action)

        # Setting the default selected electrograms
        new_system.egm_points = np.asarray([0], dtype=int)

        new_system._determine_available_fields()
        self.interpolate_openEP_fields(system=new_system)
        self.add_view(new_system)

    def _load_openCARP(self):
        """
        Load a set of openCARP files.

        Requires that both a *.pts and *.elem file are selected (and no other file).
        Auxillary added can later be by going to 'File > Add data to'.
        """

        dialogue = QtWidgets.QFileDialog()
        dialogue.setWindowTitle('Load a set of openCARP files')
        dialogue.setDirectory(QtCore.QDir.currentPath())
        dialogue.setFileMode(QtWidgets.QFileDialog.ExistingFiles)
        dialogue.setNameFilter("openCARP files (*.pts *.elem)")

        error = QtWidgets.QMessageBox()
        error.setIcon(QtWidgets.QMessageBox.Critical)
        error.setText("File selection error")
        error.setInformativeText(
            "Please select a single set of openCARP files (one each of *.pts, *.elem)."
        )
        error.setWindowTitle("Error")

        if dialogue.exec_():

            filenames = dialogue.selectedFiles()
            extensions = [pathlib.Path(f).suffix for f in filenames]

            if len(filenames) != 2:
                error.exec_()
                return
            elif (".pts" not in extensions) or (".elem" not in extensions):
                error.exec_()
                return

            self._initialise_openCARP(
                points=filenames[extensions.index(".pts")],
                indices=filenames[extensions.index(".elem")],
            )

    def _initialise_openCARP(self, points, indices):
        """
        Initialise data from an openCARP simulation.

        A new system is created. If no other systems exist, this system is made the active
        system.

        Options for adding auxillary data to the system and opening a new 3D-viewer are added
        to the 'File > Add data to' and the 'View > Add view for' menus, respectively.
        """

        carp = openep.load_openCARP(
            points=points,
            indices=indices,
        )

        # we need a unique key for each system
        while str(self._system_counter) in self.systems:
            self._system_counter += 1

        new_system = openep.view.system_manager.System(
            name=str(self._system_counter),
            folder=pathlib.Path(points).parent.resolve(),
            type="openCARP",
            data=carp,
        )

        self.systems[new_system.name] = new_system
        self._system_counter += 1
        self._active_system = new_system if self._active_system is None else self._active_system

        self.update_system_manager_table(new_system)

        # We need to dynamically add options for loading data/creating 3D viewers to the main menubar
        add_data_to_system_menu = QtWidgets.QMenu(points, self)
        add_unipolar_action = QtWidgets.QAction("Unipolar electrograms", self)
        add_unipolar_action.triggered.connect(lambda: self._add_data_to_openCARP(new_system))
        add_data_to_system_menu.addAction(add_unipolar_action)
        self.system_main.add_data_menu.addMenu(add_data_to_system_menu)

        add_view_for_system_action = QtWidgets.QAction(points, self)
        add_view_for_system_action.triggered.connect(lambda: self.add_view(new_system))
        self.system_main.add_view_menu.addAction(add_view_for_system_action)

        # Setting the default selected electrograms
        new_system.egm_points = np.asarray([0], dtype=int)

    def _add_data_to_openCARP(self, system):
        """
        Add data into an CARPData object.

        If there are no plotters for this system, a new 3D viewer will automatically be created.
        """

        dialogue = QtWidgets.QFileDialog()
        dialogue.setWindowTitle('Add an openCARP data file')
        dialogue.setDirectory(QtCore.QDir.currentPath())
        dialogue.setFileMode(QtWidgets.QFileDialog.ExistingFile)
        dialogue.setNameFilter("openCARP data file (*.dat)")

        if dialogue.exec_():

            filename = dialogue.selectedFiles()[0]
            system.data.add_data('unipolar egm', filename)
            if len(system.plotters) == 0:
                self.add_view(system)

    def _export_data_to_openCARP(self, system):
        """Export mesh data form an OpenEP dataset into openCARP format."""

        dialogue = QtWidgets.QFileDialog()
        dialogue.setWindowTitle('Save As (prefix only)')
        dialogue.setDirectory(QtCore.QDir.currentPath())
        dialogue.setFileMode(QtWidgets.QFileDialog.AnyFile)
        dialogue.setNameFilter("No extension ()")

        prefix, _ = QtWidgets.QFileDialog.getSaveFileName(
            self,
            caption="Save As (prefix only)",
            directory=QtCore.QDir.currentPath(),
            filter="Prefix only ()",
        )

        error = QtWidgets.QMessageBox()
        error.setIcon(QtWidgets.QMessageBox.Critical)
        error.setText("File selection error")
        error.setInformativeText(
            "Please select a filename prefix - do not specify the extension."
        )
        error.setWindowTitle("Error")

        prefix_path = pathlib.Path(prefix).resolve()
        if prefix_path.suffix:
            error.exec_()
            return

        openep.export_openCARP(
            case=system.data,
            prefix=prefix,
        )

    def add_view(self, system):
        """
        Create a new CustomDockWidget for the given system (i.e. open a new 3D viewer).

        By default, the bipolar voltage is drawn.
        """

        if system.data.electric is None:
            error = QtWidgets.QMessageBox()
            error.setIcon(QtWidgets.QMessageBox.Critical)
            error.setText("No Data Error")
            error.setInformativeText(
                "Please load data into the system before creating a new view.\n"
                "File > Add data to"
            )
            error.setWindowTitle("Error")
            error.exec_()
            return

        dock, plotter, mesh = system.create_dock()
        add_mesh_kws = system._create_default_kws()
        free_boundaries = openep.mesh.get_free_boundaries(mesh)

        index = len(system.plotters)
        plotter.lower_limit.returnPressed.connect(lambda: self.update_colourbar_limits(system, index=index))
        plotter.upper_limit.returnPressed.connect(lambda: self.update_colourbar_limits(system, index=index))
        plotter.opacity.valueChanged.connect(lambda: self.update_opacity(system, index=index))

        system.docks.append(dock)
        system.plotters.append(plotter)
        system.meshes.append(mesh)
        system.add_mesh_kws.append(add_mesh_kws)
        system.free_boundaries.append(free_boundaries)

        # We can't put this into a for loop.
        # For some reason, when iterating over the items in plotter.scalar_field_actions,
        # all actions are passed the key in the iteration for the scalars argument
        if 'Bipolar voltage' in plotter.scalar_field_actions:
            plotter.scalar_field_actions['Bipolar voltage'].triggered.connect(
                lambda: self.change_active_scalars(system, index=index, scalars='Bipolar voltage')
            )
        if 'Unipolar voltage' in plotter.scalar_field_actions:
            plotter.scalar_field_actions['Unipolar voltage'].triggered.connect(
                lambda: self.change_active_scalars(system, index=index, scalars='Unipolar voltage')
            )
        if 'Clinical bipolar voltage' in plotter.scalar_field_actions:
            plotter.scalar_field_actions['Clinical bipolar voltage'].triggered.connect(
                lambda: self.change_active_scalars(system, index=index, scalars='Clinical bipolar voltage')
            )
        if 'Clinical unipolar voltage' in plotter.scalar_field_actions:
            plotter.scalar_field_actions['Clinical unipolar voltage'].triggered.connect(
                lambda: self.change_active_scalars(system, index=index, scalars='Clinical unipolar voltage')
            )
        if 'Clinical LAT' in plotter.scalar_field_actions:
            plotter.scalar_field_actions['Clinical LAT'].triggered.connect(
                lambda: self.change_active_scalars(system, index=index, scalars='Clinical LAT')
            )
        if 'Clinical force' in plotter.scalar_field_actions:
            plotter.scalar_field_actions['Clinical force'].triggered.connect(
                lambda: self.change_active_scalars(system, index=index, scalars='Clinical force')
            )
        if 'Clinical impedance' in plotter.scalar_field_actions:
            plotter.scalar_field_actions['Clinical impedance'].triggered.connect(
                lambda: self.change_active_scalars(system, index=index, scalars='Clinical impedance')
            )

        dock.setWindowTitle(f"{system.name}: {mesh.active_scalars_info.name}")
        active_scalars_name = mesh.active_scalars_info.name
        self.draw_map(
            mesh=mesh,
            plotter=plotter,
            data=None,
            add_mesh_kws=add_mesh_kws,
            free_boundaries=free_boundaries,
        )

        # Doing add_mesh and setting a title to the colour bar causes pyvista to add
        # the data to the point_data array with the name of the title, and then sets this
        # new point_data array as the active scalars
        mesh.set_active_scalars(active_scalars_name)

        # If this is the first 3d viewer for the first system loaded, we need to update the egms etc.
        if (system.name == self._active_system.name) and (len(system.plotters) == 1):
            self.change_active_system(system)

        self.addDockWidget(Qt.LeftDockWidgetArea, dock)
        self.highlight_menubars_of_active_plotters()

    def update_colourbar_limits(self, system, index):
        """Update colourbar limits for a given plotter.

        Args:
            system (System): system containing the plotter whose colourbar
                limits will be changed
            index (int): index of plotter to modify
        """

        plotter = system.plotters[index]

        if not plotter.lower_limit.text().strip():
            return
        elif not plotter.upper_limit.text().strip():
            return

        lower_limit = float(system.plotters[index].lower_limit.text())
        upper_limit = float(system.plotters[index].upper_limit.text())
        system.add_mesh_kws[index]["clim"] = [lower_limit, upper_limit]

        actor = plotter.renderer._actors['mesh']
        actor_mapper = actor.GetMapper()
        actor_mapper.SetScalarRange((lower_limit, upper_limit))

    def update_opacity(self, system, index):
        """Update the opacity for a given ploter.

        Args:
            system (System): system containing the plotter whose opacity
                will be changed
            index (int): index of plotter to modify
        """

        plotter = system.plotters[index]
        system.add_mesh_kws[index]['opacity'] = plotter.opacity.value()

        actor = plotter.renderer._actors['mesh']
        actor_properties = actor.GetProperty()
        actor_properties.SetOpacity(plotter.opacity.value())

    def change_active_system(self, system):
        """
        Update selected active system.

        Make the menubars blue for all 3d viewers of the active system.
        If the active system has electrograms, plot these in the electrogram viewer.
        Change the window of interest range slider in the electrogram viewer to reflect
        those of the active system.
        """

        self._active_system = system

        self.highlight_menubars_of_active_plotters()

        self.remove_woi_slider()
        self.check_for_available_electrograms()
        if self.has_electrograms:

            self.egm_axis.axis('on')  # make sure we can see the axes now
            self.egm_axis.set_xlim(self.egm_times[0]-100, self.egm_times[-1]+100)
            self.egm_axis.set_ylim(-1, 13)

            # We can't create the slider earlier because we need to know the min and max egm times
            self.create_woi_slider()
            self.initialise_woi_slider_limits()
            self.egm_slider.on_changed(self.update_woi_slider_limits)
            self.update_electrograms_from_stored()

        else:
            self.egm_axis.cla()
            self.egm_axis.axis('off')
            self.egm_canvas.draw()

    def change_active_scalars(self, system, index, scalars):
        """Cahnge the scalar values that are being projected onto the mesh."""

        dock = system.docks[index]
        dock.setWindowTitle(f"{system.name}: {scalars}")

        mesh = system.meshes[index]
        mesh.set_active_scalars(name=scalars)

    def update_scalar_fields(self, event):
        """
        If the active system is an OpenEP case:
        Interpolate EGM data onto the surface and draw a map if necessary.

        Or, if the active system is an openCARP simulation:
        Recalculate the bipolar and unipolar voltages using the new window of interest.

        The event argument is ignored. It is there because when
        self.set_woi_button.on_clicked (mpl.widgets.Button) is pressed,
        matplotlib passes an event to the called function (i.e. this one).
        """

        system = self._active_system
        if system.type == "OpenEP":
            self.interpolate_openEP_fields()
        elif system.type == "openCARP":
            self.recalculate_openCARP_fields()

    def interpolate_openEP_fields(self, system=None):
        """
        Interpolate EGM data onto the surface.

        We use a buffer of zero as the window of interest if determined by the user,
        via self.slider (mpl.widgets.RangeSlider) and self.set_woi_button (mpl.widgets.Button).
        """

        system = self._active_system if system is None else system
        case = system.data

        try:
            bipolar_voltage = openep.case.interpolate_voltage_onto_surface(
                case,
                max_distance=None,
                buffer=0,
            )
            unipolar_voltage = openep.case.interpolate_voltage_onto_surface(
                case,
                max_distance=None,
                buffer=0,
                bipolar=False,
            )
        except ValueError as e:

            if str(e) != "cannot reshape array of size 0 into shape (0,newaxis)":
                raise

            error = QtWidgets.QMessageBox()
            error.setIcon(QtWidgets.QMessageBox.Critical)
            error.setText("Range Error")
            error.setInformativeText(
                "There are no mapping points within the selected window of interest.\n"
                "Voltages from mapping points cannot be interpolated onto the surface."
            )
            error.setWindowTitle("Error")
            error.exec_()
            return

        # We need to modify these value so docks/plotters/meshes create in the future
        # use the correct (current) values
        system.scalar_fields['Bipolar voltage'][:] = bipolar_voltage
        system.scalar_fields['Unipolar voltage'][:] = unipolar_voltage

        # We index into the arrays so we modify the view without redrawing the scene
        for mesh in system.meshes:
            mesh.point_data['Bipolar voltage'][:] = bipolar_voltage
            mesh.point_data['Unipolar voltage'][:] = unipolar_voltage

    def recalculate_openCARP_fields(self):
        """Recalculate the scalar values that are being projected onto the mesh for an openCARP dataset."""

        system = self._active_system
        carp = system.data

        # TODO: If bipolar voltages are calculated from unipolar electrograms, the bipolar electrograms
        #       should first be reproducted using carp.bipolar_from_unipolar(electrograms[:, woi])

        if self.egm_unipolar_A_checkbox.isEnabled():

            unipolar_voltage = openep.case.calculate_voltage_from_electrograms(
                carp, buffer=0, bipolar=False,
            )

            carp.electric.unipolar_egm.voltage[:] = unipolar_voltage
            system.scalar_fields['Unipolar voltage'][:] = unipolar_voltage[:, 0]
            for mesh in system.meshes:
                mesh.point_data['Unipolar voltage'][:] = unipolar_voltage[:, 0]

        if self.egm_bipolar_checkbox.isEnabled():

            bipolar_voltage = openep.case.calculate_voltage_from_electrograms(
                carp, buffer=0, bipolar=True,
            )

            carp.electric.bipolar_egm.voltage[:] = bipolar_voltage
            system.scalar_fields['Bipolar voltage'][:] = bipolar_voltage
            for mesh in system.meshes:
                mesh.point_data['Bipolar voltage'][:] = bipolar_voltage

    def draw_map(self, mesh, plotter, data, add_mesh_kws, free_boundaries=None):
        """
        Project scalar values onto the surface of the mesh.

        Also draw the free boundaries (anatomical structures) as black lines.
        """

        plotter = openep.draw.draw_map(
            mesh=mesh,
            field=data,
            plotter=plotter,
            free_boundaries=False,
            add_mesh_kws=add_mesh_kws,
        )

        if free_boundaries is not None:
            plotter = openep.draw.draw_free_boundaries(
                free_boundaries,
                colour="black",
                width=5,
                plotter=plotter,
            )

    def check_for_available_electrograms(self):
        """Update the electrom types that can be plotted."""

        if self._active_system.data.electric is None:
            self.has_electrograms = False
            self.egm_dock.setEnabled(False)
            self.egm_reference_checkbox.setEnabled(False)
            self.egm_bipolar_checkbox.setEnabled(False)
            self.egm_unipolar_A_checkbox.setEnabled(False)
            self.egm_unipolar_B_checkbox.setEnabled(False)
            return

        electric = self._active_system.data.electric

        self.has_electrograms = False
        if electric.reference_egm is not None and electric.reference_egm.egm is not None:
            self.has_electrograms = True
            self.egm_reference_checkbox.setEnabled(True)
            self.egm_times = np.arange(electric.reference_egm.egm.shape[1])
        else:
            self.egm_reference_checkbox.setEnabled(False)

        if electric.bipolar_egm is not None and electric.bipolar_egm.egm is not None:
            self.has_electrograms = True
            self.egm_bipolar_checkbox.setEnabled(True)
            self.egm_times = np.arange(electric.bipolar_egm.egm.shape[1])
        else:
            self.egm_bipolar_checkbox.setEnabled(False)

        if electric.unipolar_egm is not None and electric.unipolar_egm.egm is not None:
            self.has_electrograms = True
            self.egm_unipolar_A_checkbox.setEnabled(True)
            self.egm_unipolar_B_checkbox.setEnabled(True)
            self.egm_times = np.arange(electric.unipolar_egm.egm.shape[1])
        else:
            self.egm_unipolar_A_checkbox.setEnabled(False)
            self.egm_unipolar_B_checkbox.setEnabled(False)

        if self.has_electrograms:
            self.egm_dock.setEnabled(True)
        else:
            self.egm_dock.setEnabled(False)

    def update_electrograms_from_text(self):
        """
        Extract electrograms at specified indices and re-plot.

        The points are specified by the user and set via use of the
        self.egm_select (QLineEdit) and button_egm_select (QPushButton) widgets.
        """

        # Get data for new set of points
        egm_text = self.egm_selection.text()
        selected_indices = self.egm_selection_filter.findall(egm_text)
        self._active_system.egm_points = np.asarray(selected_indices, dtype=int)

        self.extract_electrograms()
        self.plot_electrograms()

    def update_electrograms_from_stored(self):

        # Retrieve data for the stored set of points
        egm_points = self._active_system.egm_points
        egm_text = " ".join(egm_points.astype(str))
        self.egm_selection.setText(egm_text)

        self.extract_electrograms()
        self.plot_electrograms()

    def extract_electrograms(self):

        system = self._active_system

        if self.egm_reference_checkbox.isEnabled():
            self.egm_reference_traces, self.egm_names = openep.case.get_electrograms_at_points(
                system.data,
                within_woi=False,
                buffer=0,
                indices=system.egm_points,
                egm_type="reference",
                return_names=True,
                return_lat=False,
            )

        if self.egm_bipolar_checkbox.isEnabled():
            self.egm_bipolar_traces, self.egm_names = openep.case.get_electrograms_at_points(
                system.data,
                within_woi=False,
                buffer=0,
                indices=system.egm_points,
                egm_type="bipolar",
                return_names=True,
                return_lat=False,
            )

        if self.egm_unipolar_A_checkbox.isEnabled():
            unipolar_traces, self.egm_names = openep.case.get_electrograms_at_points(
                system.data,
                within_woi=False,
                buffer=0,
                indices=system.egm_points,
                egm_type="unipolar",
                return_names=True,
                return_lat=False,
            )
            self.egm_unipolar_A_traces = unipolar_traces[:, :, 0]
            self.egm_unipolar_B_traces = unipolar_traces[:, :, 1]

    def plot_electrograms(self):
        """
        Plot electrograms for the currently-selected set of points.

        Electrograms must first have been extracted using update_electrograms.
        Here we will plot the reference, bipolar, unipolar A, unipolar B
        electrograms for each point.
        """

        # Set up axis for new plots
        ylim = self.egm_axis.get_ylim()
        xlim = self.egm_axis.get_xlim()

        self.egm_axis.cla()
        self.egm_axis.set_yticklabels([])
        self.egm_axis.set_ylim(ylim)
        self.egm_axis.set_xlim(xlim)

        # Reference voltage
        if self.egm_reference_checkbox.isEnabled() and self.egm_reference_checkbox.isChecked():

            _, self.egm_axis.axes = openep.draw.plot_electrograms(
                self.egm_times,
                self.egm_reference_traces,
                axis=self.egm_axis.axes,
                colour="xkcd:scarlet",
                y_separation=2,
            )

        # Bipolar voltage
        if self.egm_bipolar_checkbox.isEnabled() and self.egm_bipolar_checkbox.isChecked():

            _, self.egm_axis.axes = openep.draw.plot_electrograms(
                self.egm_times,
                self.egm_bipolar_traces,
                axis=self.egm_axis.axes,
                colour="xkcd:cerulean",
                y_separation=2,
            )

        # Unipolar A voltage
        if self.egm_unipolar_A_checkbox.isEnabled() and self.egm_unipolar_A_checkbox.isChecked():

            _, self.egm_axis.axes = openep.draw.plot_electrograms(
                self.egm_times,
                self.egm_unipolar_A_traces,
                axis=self.egm_axis.axes,
                colour="xkcd:tree green",
                y_start=1,
                y_separation=2,
            )

        # Unipolar B voltage
        if self.egm_unipolar_B_checkbox.isEnabled() and self.egm_unipolar_B_checkbox.isChecked():

            _, self.egm_axis.axes = openep.draw.plot_electrograms(
                self.egm_times,
                self.egm_unipolar_B_traces,
                axis=self.egm_axis.axes,
                colour="xkcd:pumpkin",
                y_start=1,
                y_separation=2,
            )

        # Add labels if necessary
        yticks = []
        yticklabels = []
        separations = np.arange(self._active_system.egm_points.size) * 2

        if (
            (self.egm_reference_checkbox.isEnabled() and self.egm_reference_checkbox.isChecked()) or
            (self.egm_bipolar_checkbox.isEnabled() and self.egm_bipolar_checkbox.isChecked())
        ):
            yticks.extend(separations)
            yticklabels.extend(self.egm_names)

        # Unipolar A and B are shifted above the bipolar and reference electrograms for clarity
        if (
            (self.egm_unipolar_A_checkbox.isEnabled() and self.egm_unipolar_A_checkbox.isChecked()) or
            (self.egm_unipolar_B_checkbox.isEnabled() and self.egm_unipolar_B_checkbox.isChecked())
        ):
            yticks.extend(separations + 1)
            yticklabels.extend(self.egm_names)

        self.egm_axis.set_yticks(yticks)
        self.egm_axis.set_yticklabels(yticklabels)

        # draw vertical lines at the window of interest
        self.initialise_woi_slider_limits()

        self.egm_canvas.draw()

    def create_woi_slider(self):
        """Add a window of interest range slider to the electrogram canvas"""

        self.egm_slider = openep.view.canvases.add_woi_range_slider(
            figure=self.egm_figure,
            axis=[0.1535, 0.05, 0.7185, 0.01],  # location of the RangeSlider on the figure
            valmin=self.egm_times[0],
            valmax=self.egm_times[-1],
        )

    def remove_woi_slider(self):
        """Remove the window of interest range slider from the electrogram canvas"""

        if self.egm_slider is not None:
            self.egm_slider.ax.remove()
            self.egm_slider = None

    def initialise_woi_slider_limits(self):
        """Set the limits to be the window of interest."""

        annotations = self._active_system.data.electric.annotations

        start_woi, stop_woi = annotations.window_of_interest[0]
        start_woi += annotations.reference_activation_time[0]
        stop_woi += annotations.reference_activation_time[0]

        self.egm_slider.set_val([start_woi, stop_woi])

        self.egm_slider_lower_limit = self.egm_axis.axvline(
            start_woi,
            color="grey",
            linestyle='--',
            linewidth=0.8,
            alpha=0.6,
        )

        self.egm_slider_upper_limit = self.egm_axis.axvline(
            stop_woi,
            color="grey",
            linestyle='--',
            linewidth=0.8,
            alpha=0.6,
        )

    def update_woi_slider_limits(self, val):
        """
        Take the min and max values from the RangeSlider widget.
        Use this to set the window of interest and to change the x location
        of the two axvlines drawn on the EGM canvas.
        """

        annotations = self._active_system.data.electric.annotations

        # from https://matplotlib.org/devdocs/gallery/widgets/range_slider.html
        start_woi, stop_woi = val
        self.egm_slider_lower_limit.set_xdata([start_woi, start_woi])
        self.egm_slider_upper_limit.set_xdata([stop_woi, stop_woi])

        reference_annotation = annotations.reference_activation_time[0]
        annotations.window_of_interest[:, 0] = start_woi - reference_annotation
        annotations.window_of_interest[:, 1] = stop_woi - reference_annotation

    def highlight_menubars_of_active_plotters(self):
        """Make the menubars of all docks in the active system blue."""

        for system_name, system in self.systems.items():
            if system_name == self._active_system.name:
                for dock in system.docks:
                    dock.main.menubar.setStyleSheet(
                        "QMenuBar{"
                        "background-color: #5a86ad;"  # xkcd:dusty blue
                        "color: black;"
                        "border: None;"
                        "}"
                    )
                # The below can be used to put a red border around the plotter (but not the entire window)
                # for plotter in system.plotters:
                #     plotter.renderer.add_border(color="red", width=5)
            else:
                # The below can be used to remove the red border around the plotter
                # for plotter in system.plotters:
                #     plotter.renderer.add_border(color=None, width=5)
                for dock in system.docks:
                    dock.main.menubar.setStyleSheet(
                        "QMenuBar{"
                        "background-color: #95a3a6;"  # xkcd:cool grey
                        "color: black;"
                        "border: None;"
                        "}"
                    )


def main():

    # Create an instance of Qapplication
    app = QtWidgets.QApplication(sys.argv)
    app.setWindowIcon(QtGui.QIcon(openep.view.images.LOGO))

    # Create an instance of GUI
    window = OpenEPGUI()
    window.showMaximized()

    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
