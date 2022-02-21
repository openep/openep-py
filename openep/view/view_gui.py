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
from distutils.version import LooseVersion
import platform
import sys
import os
import pathlib

import qdarkstyle
from PySide2 import QtCore, QtGui, QtWidgets
from PySide2.QtCore import Qt
from PySide2.QtWidgets import QAction
from matplotlib.backend_bases import MouseButton
import numpy as np

import openep

import openep.view.custom_widgets
import openep.view.mapping_points
import openep.view.annotations_ui
import openep.view.plotters_ui
import openep.view.system_manager_ui
import openep.view.system_ui

import openep.view.plotters
import openep.view.system_manager

import openep.view.static


# Fixes issues with window not displaying on Big Sur
# https://bugreports.qt.io/browse/QTBUG-87014
if (sys.platform == 'darwin' and
        LooseVersion(platform.mac_ver()[0]) >= LooseVersion("10.16") and
        LooseVersion(QtCore.qVersion()) < LooseVersion("6") and
        "QT_MAC_WANTS_LAYER" not in os.environ):
            os.environ["QT_MAC_WANTS_LAYER"] = "1"


class OpenEPGUI(QtWidgets.QMainWindow):

    def __init__(self):

        super().__init__()

        self._init_systems()
        self._init_ui()
        self._create_system_manager_ui()
        self._create_annotate_dock()
        self._create_mapping_points_dock()
        self._add_dock_widgets()
        self._disable_dock_widgets()

    def _init_systems(self):
        """Containers and variables for keeping track of the systems loaded into the GUI."""

        self.system_manager = openep.view.system_manager.SystemManager()
        self.egm_times = None

    def _init_ui(self):
        """
        Initialise the GUI.

        Lots of the layout is handled automatically by QtWidgets.QMainWindow.
        """

        self.setMinimumSize(800, 600)
        self.setWindowTitle("OpenEP: The open-source solution for electrophysiology data analysis")
        self.setToolTipDuration(-1)
        self.setAttribute(QtCore.Qt.WA_MacOpaqueSizeGrip, True)

    def _create_system_manager_ui(self):
        """Create a dockable widget that for managing systems loaded into the GUI"""

        self.system_manager_ui = openep.view.system_manager_ui.SystemManagerDockWidget(
            title="Systems"
        )

        # We also need to connect up the actions
        self.system_manager_ui.main.load_openep_mat_action.triggered.connect(self._load_openep_mat)
        self.system_manager_ui.main.load_opencarp_action.triggered.connect(self._load_opencarp)

    def _create_annotate_dock(self):
        """
        Create a dockable widget for annotating electrograms with matplotlib.

        The user can use a dropdown menu to select which electrogram to annotate.
        All electrical data for that point will be plotted - electrograms, activation
        time etc.
        """
        
        self.annotate_dock = openep.view.annotations_ui.AnnotationWidget(title="Annotate")
    
        self.annotate_dock.egm_selection.currentIndexChanged[int].connect(self.update_annotation_plot)
        self.annotate_dock.canvas.mpl_connect('key_press_event', self.annotation_on_key_press)
        self.annotate_dock.canvas.mpl_connect('scroll_event', self.annotation_on_scroll_wheel)  # this is scrolling with the wheel, not scrollbar
        
        # We will need to disconnect and reconnect these events through user interaction
        self.annotate_dock.cid_button_press_event = self.annotate_dock.canvas.mpl_connect(
            'button_press_event',
            self.annotation_on_button_press,
        )
        self.annotate_dock.cid_motion_notify_event_cursor_style = self.annotate_dock.canvas.mpl_connect(
            'motion_notify_event',
            self.annotate_dock._on_mouse_move_cursor_style,
        )
        
        # remove all default key bindings
        # see https://stackoverflow.com/a/35631062/17623640
        # self.annotate_dock.canvas.mpl_disconnect(self.annotate_dock.canvas.manager.key_press_handler_id)
        # But this doesn't work for some reason
        
        # hide the figure until we have data to plot
        self.annotate_dock.deactivate_figure()
    
    def _create_mapping_points_dock(self):
        """Create a dockable widget for displaying mapping points info."""

        self.mapping_points = openep.view.mapping_points.MappingPointsDock(
            title="Mapping points",
            system=None,
        )
        self.mapping_points.delete_points.triggered.connect(
            lambda: self.delete_mapping_points(table=self.mapping_points.table, restore=False)
        )
        self.mapping_points.table.horizontalHeader().sortIndicatorChanged.connect(
            lambda: self.sort_mapping_points_and_recycle_bin(self.mapping_points, self.recycle_bin),
        )
        self.mapping_points._ignore_sort_signal = False

        self.recycle_bin = openep.view.mapping_points.RecycleBinDock(
            title="Recycle bin",
            system=None,
            model=self.mapping_points.model,
            proxy_model=self.mapping_points.proxy_model,
        )
        self.recycle_bin.restore_points.triggered.connect(
            lambda: self.delete_mapping_points(table=self.recycle_bin.table, restore=True)
        )

        self.recycle_bin.table.horizontalHeader().sortIndicatorChanged.connect(
            lambda: self.sort_mapping_points_and_recycle_bin(self.recycle_bin, self.mapping_points)
        )
        self.recycle_bin._ignore_sort_signal = False

    def _add_dock_widgets(self):
        """
        Add dockable widgets to the main window.
        """

        self.addDockWidget(Qt.RightDockWidgetArea, self.system_manager_ui)
        self.splitDockWidget(self.system_manager_ui, self.mapping_points, Qt.Vertical)
        self.splitDockWidget(self.mapping_points, self.annotate_dock, Qt.Vertical)
        self.splitDockWidget(self.mapping_points, self.recycle_bin, Qt.Horizontal)
        
        for dock in [self.system_manager_ui, self.annotate_dock, self.mapping_points, self.recycle_bin]:
            dock.setAllowedAreas(Qt.AllDockWidgetAreas)

        self.setDockOptions(self.GroupedDragging | self.AllowTabbedDocks | self.AllowNestedDocks)
        self.setTabPosition(Qt.AllDockWidgetAreas, QtWidgets.QTabWidget.North)

    def _disable_dock_widgets(self):
        """
        Ignore all key presses in the dock widgets.

        This is required if there is no file loaded, otherwise
        the GUI will crash.
        """

        self.annotate_dock.setEnabled(False)

    def _load_openep_mat(self):
        """
        Load an OpenEP case.

        Currently, only MATLAB files are supported.
        """

        dialogue = QtWidgets.QFileDialog()
        dialogue.DontUseNativeDialog = True
        dialogue.setWindowTitle('Load an OpenEP file')
        dialogue.setDirectory(QtCore.QDir.currentPath())
        dialogue.setFileMode(QtWidgets.QFileDialog.ExistingFile)
        dialogue.setNameFilter("MATLAB file (*.mat)")

        if dialogue.exec_():
            filename = dialogue.selectedFiles()[0]
            self._initialise_openep_mat(filename)

    def _initialise_openep_mat(self, filename):
        """Initialise data from an OpenEP case object.

        Create separate meshes for the two BackgroundPlotters.
        Interpolate data ready to be plotted.
        """

        case = openep.load_openep_mat(filename)
        basename = str(pathlib.Path(filename).resolve().with_suffix(""))
        system = self.system_manager.add_system(
            case=case,
            basename=basename,
            type="OpenEP",
        )
        self.update_system_manager_table(system)

        # We need to dynamically add options for exporting data/creating 3D viewers to the main menubar
        export_action = self.system_manager_ui.create_export_action(
            system_basename=system.basename,
            export_name="as openCARP",
        )
        view_action = self.system_manager_ui.create_view_action(
            system_basename=system.basename,
        )

        export_action.triggered.connect(lambda: self._export_data_to_openCARP(system))
        view_action.triggered.connect(lambda: self.add_view(system))

        self.interpolate_openEP_fields(system=system)
        self.add_view(system)

    def _load_opencarp(self):
        """
        Load a set of openCARP files.

        Requires that both a *.pts and *.elem file are selected (and no other file).
        Auxillary added can later be by going to 'File > Add data to'.
        """

        dialogue = QtWidgets.QFileDialog()
        dialogue.DontUseNativeDialog = True
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

        case = openep.load_opencarp(
            points=points,
            indices=indices,
        )
        basename = str(pathlib.Path(points).resolve().with_suffix(""))
        system = self.system_manager.add_system(
            case=case,
            basename=basename,
            type="openCARP",
        )
        self.update_system_manager_table(system)

        # We need to dynamically add option for loading data to the main menubar.
        # But we shouldn't an option to create 3d viewer yet because we have no scalar fields at the minute.
        # Instead, this is done in self._add_data_to_openCARP
        add_data_action = self.system_manager_ui.create_add_data_action(
            system_basename=system.basename,
            data_type="Unipolar electrograms",
        )
        add_data_action.triggered.connect(lambda: self._add_data_to_openCARP(system=system, data_type='unipolar_egm'))
        # TODO: support loading of different data to the system, e.g. arbitrary scalar fields that can be named by the user

    def _add_data_to_openCARP(self, system, data_type):
        """
        Add data into an CARPData object.

        If there are no plotters for this system, a new 3D viewer will automatically be created.
        """

        dialogue = QtWidgets.QFileDialog()
        dialogue.DontUseNativeDialog = True
        dialogue.setWindowTitle('Add an openCARP data file')
        dialogue.setDirectory(QtCore.QDir.currentPath())
        dialogue.setFileMode(QtWidgets.QFileDialog.ExistingFile)
        dialogue.setNameFilter("openCARP data file (*.dat)")

        if dialogue.exec_():

            filename = dialogue.selectedFiles()[0]

            if data_type == "unipolar_egm":
                unipolar = np.loadtxt(filename)
                system.case.add_unipolar_electrograms(
                    unipolar=unipolar,
                    add_bipolar=True,
                    add_annotations=True,
                )

            system._determine_available_fields()  # ensure we can access the loaded data

            if len(system.plotters) == 0:

                # We need to dynamically add an option for creating 3D viewers to the main menubar
                add_view_for_system_action = QAction(str(system.basename), self)
                add_view_for_system_action.triggered.connect(lambda: self.add_view(system))
                self.system_manager_ui.main.add_view_menu.addAction(add_view_for_system_action)

                # And also create the first 3d viewer
                self.add_view(system)

    def _export_data_to_openCARP(self, system):
        """Export mesh data form an OpenEP dataset into openCARP format."""

        dialogue = QtWidgets.QFileDialog()
        dialogue.DontUseNativeDialog = True
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
            case=system.case,
            prefix=prefix,
        )

    def update_system_manager_table(self, system):
        """Update the System manager table when a new system is loaded."""

        row_number = len(self.system_manager.systems)

        system.name_widget, *_, active_widget = self.system_manager_ui.update_system_manager_table(
            name=str(system.name),
            basename=str(system.basename),
            data_type=system.type,
            is_active=True if self.system_manager.active_system.name == system.name else False,
            row_number=row_number,
        )

        # We need to connect up the actions
        system.name_widget.returnPressed.connect(lambda: self.update_system_name(system))
        active_widget.clicked.connect(lambda: self.change_active_system(system))

    def update_system_name(self, system):
        """
        Change the name of a system.

        We also need to change the window title of its dock widgets.
        """

        new_name = system.name_widget.text()

        try:
            self.system_manager.update_system_name(system=system, new_name=new_name)

        except KeyError:

            system.name_widget.setText(system.name)

            error = QtWidgets.QMessageBox()
            error.setIcon(QtWidgets.QMessageBox.Critical)
            error.setText("Name Error")
            error.setInformativeText(
                "System names must be unique."
            )
            error.setWindowTitle("Error")
            error.exec_()

        else:
            for dock, mesh in zip(system.docks, system.meshes):
                dock.setWindowTitle(f"{system.name}: {mesh.active_scalars_info.name}")

    def add_view(self, system):
        """Create a new CustomDockWidget for the given system (i.e. open a new 3D viewer)."""

        if system.case.electric.unipolar_egm.voltage is None:
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

        plotter = openep.view.plotters.create_plotter()

        index = len(system.plotters)
        plotter_is_secondary_view = index > 0
        plotter_widgets = openep.view.plotters_ui.create_plotter_layout(
            plotter=plotter,
            add_link_views_button=plotter_is_secondary_view,
        )
        central_widget, plotter.lower_limit, plotter.upper_limit, plotter.opacity, plotter.link_view_with_primary = \
            plotter_widgets

        dock = openep.view.system_ui.create_system_dock(central_widget=central_widget)

        # Add a Field menu to the menubar. This is used for selecting the scalar field to project onto the surface.
        dock, plotter = openep.view.system_ui.add_field_menu(
            dock=dock,
            plotter=plotter,
            system_name=system.name,
            scalar_fields=system.scalar_fields,
        )
        for action_name, action in plotter.scalar_field_actions.items():
            action.toggled.connect(
                lambda checked, scalars=action_name: self.change_active_scalars(
                    system, index=index, scalars=scalars
                )
            )

        # Add a Show/Hide menu to the menubar. This is used for showing/hiding the surface mesh, mapping points,
        # and surface-projected mapping points
        dock, plotter = openep.view.system_ui.add_show_menu(
            dock=dock,
            plotter=plotter,
        )
        for action_name, action in plotter.show_actions.items():
            action.toggled.connect(
                lambda checked, actor_name=action_name: self.update_actor_visibility(
                    plotter, actor_name=actor_name,
                )
            )
        
        # Connect slots for controlling colourbar limits, mesh opacity, and linking of plotters
        plotter.lower_limit.returnPressed.connect(lambda: self.update_colourbar_limits(system, index=index))
        plotter.upper_limit.returnPressed.connect(lambda: self.update_colourbar_limits(system, index=index))
        plotter.opacity.valueChanged.connect(lambda: self.update_opacity(system, index=index))
        if plotter_is_secondary_view:
            system.plotters[0].link_views_across_plotters(plotter)
            plotter.link_view_with_primary.toggled.connect(lambda: self.link_views_across_plotters(system, index=index))

        mesh = system.create_mesh()
        mapping_points = system.create_mapping_points_mesh()
        projected_discs = system.create_surface_discs_mesh()
        add_mesh_kws, add_points_kws, add_discs_kws = system._create_default_kws()
        free_boundaries = openep.mesh.get_free_boundaries(mesh)

        system.docks.append(dock)
        system.plotters.append(plotter)
        system.meshes.append(mesh)
        system.mapping_points_meshes.append(mapping_points)
        system.surface_projected_discs_meshes.append(projected_discs)
        system.add_mesh_kws.append(add_mesh_kws)
        system.free_boundaries.append(free_boundaries)

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
        # new point_data array as the active scalars. We need to undo this unwanted behaviour.
        mesh.set_active_scalars(active_scalars_name)

        self.draw_mapping_points(
            mesh=mapping_points,
            plotter=plotter,
            add_points_kws=add_points_kws,
        )
        plotter.renderer._actors['Mapping points'].SetVisibility(False)

        self.draw_nearest_points_discs(
            mesh=projected_discs,
            plotter=plotter,
            add_discs_kws=add_discs_kws,
        )
        plotter.renderer._actors['Surface-projected mapping points'].SetVisibility(False)

        # If this is the first 3d viewer for the first system loaded, we need to update the egms etc.
        if (system.name == self.system_manager.active_system.name) and (len(system.plotters) == 1):
            self.change_active_system(system)
            self.addDockWidget(Qt.LeftDockWidgetArea, dock)
        else:
            active_parent = self.system_manager.active_system.docks[0]
            self.tabifyDockWidget(active_parent, dock)
            dock.setTitleBarWidget(QtWidgets.QWidget())

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

        actor = plotter.renderer._actors['Surface']
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

        for actor_name, actor in plotter.renderer._actors.items():
            if actor_name == "Surface" or actor_name.startswith('free_boundary_'):
                actor_properties = actor.GetProperty()
                actor_properties.SetOpacity(plotter.opacity.value())

    def update_actor_visibility(self, plotter, actor_name):
        """Make an actor (e.g. a mesh or point cloud) visible or invisible."""

        show = True if plotter.show_actions[actor_name].isChecked() else False
        plotter.renderer._actors[actor_name].SetVisibility(show)

        if actor_name == "Surface":
            for actor_name, actor in plotter.renderer._actors.items():
                if actor_name.startswith('free_boundary_'):
                    actor.SetVisibility(show)

    def link_views_across_plotters(self, system, index):
        """Link or unlink the views of two plotters.

        Args:
            system (System): system containing the plotter whose view will be
                linked/unliked with the primary plotter of the system
            index (int): index of plotter to link/unlink the view of
        """

        plotter = system.plotters[index]

        if plotter.link_view_with_primary.isChecked():

            primary_plotter = system.plotters[0]
            primary_plotter.link_views_across_plotters(other_plotter=plotter)

            plotter.link_view_with_primary.setToolTip("3D viewer is linked with the primary 3D viewer")
            plotter.link_view_with_primary.setIcon(QtGui.QIcon(openep.view.static.LINK_ICON))

            return

        # TODO: submit a PR to pyvista to add a 'keep_camera_position' argument to 'plotter.unlink_views'
        original_camera_position = plotter.camera_position.to_list()
        plotter.unlink_views()
        plotter.camera_position = original_camera_position

        plotter.link_view_with_primary.setToolTip("3D viewer is independent of the primary 3D viewer")
        plotter.link_view_with_primary.setIcon(QtGui.QIcon(openep.view.static.UNLINK_ICON))

    def change_active_system(self, system):
        """
        Update selected active system.

        Make the menubars blue for all 3d viewers of the active system.
        If the active system has electrograms, plot these in the electrogram viewer.
        Change the window of interest range slider in the electrogram viewer to reflect
        those of the active system.
        """

        self.system_manager.active_system = system

        self.highlight_menubars_of_active_plotters()
        self.check_available_electrograms()
        if self.has_electrograms:
            xmin = self.egm_times[0] - 100
            xmax = self.egm_times[-1] + 100
            self.annotate_dock.activate_figure(xmin, xmax)
            self.create_annotation_plot()
            self.initialise_annotation_egm_selection()

        else:
            self.annotate_dock.deactivate_figure()
        
        self.mapping_points.model.system = system  # the recycle bin shares the same model
        self._hide_mapping_points()

    def change_active_scalars(self, system, index, scalars):
        """Change the scalar values that are being projected onto the mesh."""

        dock = system.docks[index]
        dock.setWindowTitle(f"{system.name}: {scalars}")

        mesh = system.meshes[index]
        mesh.set_active_scalars(name=scalars)

    def update_scalar_fields(self, event=None):
        """
        If the active system is an OpenEP case:
        Interpolate EGM data onto the surface and draw a map if necessary.

        Or, if the active system is an openCARP simulation:
        Recalculate the bipolar and unipolar voltages using the new window of interest.

        The event argument is ignored. It is there because when
        self.set_woi_button.on_clicked (mpl.widgets.Button) is pressed,
        matplotlib passes an event to the called function (i.e. this one).
        """

        system = self.system_manager.active_system
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

        system = self.system_manager.active_system if system is None else system
        case = system.case

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

        system = self.system_manager.active_system
        carp = system.case

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

        Args:
            mesh (pyvista.PolyData): mesh to be added to the plotter
            plotter (BackgroundPlotter): plotter to which the mesh will be added
            data (np.ndarray): scalar values that will be used to colour the triangles of the mesh
            add_mesh_kws (dict): keyword arguments to pass to `plotter.add_mesh`
            free_boundaries (FreeBoundaries): the free boundaries of the mesh to add to the plotter
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

    def draw_mapping_points(self, mesh, plotter, add_points_kws):
        """Render mapping points as 3D spheres.

        Args:
            mesh (pyvista.PolyData): mesh to be added to the plotter
            plotter (BackgroundPlotter): plotter to which the mesh will be added
            add_points_kws (dict): keyword arguments to pass to `plotter.add_mesh`
        """

        plotter.add_mesh(mesh, **add_points_kws)

    def draw_nearest_points_discs(self, mesh, plotter, add_discs_kws):
        """Render the surface-projected mapping points as 2D discs.

        Args:
            mesh (pyvista.PolyData): mesh to be added to the plotter
            plotter (BackgroundPlotter): plotter to which the mesh will be added
            add_discs_kws (dict): keyword arguments to pass to `plotter.add_mesh`
        """

        plotter.add_mesh(mesh, **add_discs_kws)

    def check_available_electrograms(self):
        """Check whether the active case has electrograms.
        
        If so, enable to annotation viewer, and determine times at which each
        sample was taken."""

        has_reference, has_bipolar, has_unipolar = self.system_manager.check_available_electrograms()
        self.has_electrograms = any([has_reference, has_bipolar, has_unipolar])
        self.annotate_dock.setEnabled(self.has_electrograms)
        self._set_electrogram_times()

    def _set_electrogram_times(self):
        """Generate the electrogram times based on the number of points in the timeseries.

        Warning
        -------
        This assumes a 1 ms period between each point in the timeseries.
        """

        self.egm_times = self.system_manager.electrogram_times()

    def initialise_annotation_egm_selection(self):
        """Set the available electrogram for annotating"""

        electric = self.system_manager.active_system.case.electric
        
        if electric.internal_names is not None:
            names = electric.internal_names
        else:
            names = np.asarray([f"P{index}" for index in range(len(electric.bipolar_egm.egm))], dtype=str)

        self.annotate_dock.initialise_egm_selection(names)

    def create_annotation_plot(self):
        """Update the data plotted in the annotation viewer."""
        
        # TODO: this will currently fail unless there is both reference
        # and bipolar electrograms as well as ecgs. Annotations
        # are also required. And the gain of each signal.
        
        # Plot electrical data for this point
        current_index = self.annotate_dock.egm_selection.currentIndex()
        electric = self.system_manager.active_system.case.electric
        
        reference = electric.reference_egm.egm[current_index]
        bipolar = electric.bipolar_egm.egm[current_index]
        ecg = electric.ecg.ecg[current_index]
        signals = np.asarray([reference, bipolar, ecg])
        
        reference_gain = electric.reference_egm.gain[current_index]
        bipolar_gain = electric.bipolar_egm.gain[current_index]
        ecg_gain = electric.ecg.gain[current_index]
        signal_gains = np.asarray([reference_gain, bipolar_gain, ecg_gain])
        
        times = self.system_manager.electrogram_times()
        labels = np.asarray(["Ref", "Bipolar", "ECG"])

        self.annotate_dock.plot_signals(times, signals, labels, signal_gains)
        
        # add the annotations
        self.annotate_dock._initialise_annotations()
        annotations = electric.annotations
        times = self.egm_times

        local_annotation = annotations.local_activation_time[current_index]
        local_annotation_index = np.searchsorted(times, local_annotation)
        local_annotation_index = min(local_annotation_index, times.size)
        self.annotate_dock.update_annotation(
            signal=self.annotate_dock.signal_artists['Bipolar'],
            annotation=self.annotate_dock.annotation_artists['local_annotation_point'],
            annotation_line=self.annotate_dock.annotation_artists['local_annotation_line'],
            index=local_annotation_index,
        )

        reference_annotation = annotations.reference_activation_time[current_index]
        reference_annotation_index = np.searchsorted(times, reference_annotation)
        reference_annotation_index = min(reference_annotation_index, times.size)
        self.annotate_dock.update_annotation(
            signal=self.annotate_dock.signal_artists['Ref'],
            annotation=self.annotate_dock.annotation_artists['reference_annotation_point'],
            annotation_line=self.annotate_dock.annotation_artists['reference_annotation_line'],
            index=reference_annotation_index,
        )

        start_woi, stop_woi = annotations.window_of_interest[current_index]
        start_woi += reference_annotation
        stop_woi += reference_annotation
        self.annotate_dock.update_window_of_interest(start_woi, stop_woi)

        self.annotate_dock._initialise_scrollbar(start_woi=start_woi)
        self.annotate_dock.canvas.draw()
        self.annotate_dock.update_active_artist()  # ensure the annotations/woi are added to the background

    def update_annotation_plot(self):
                
        # TODO: this will currently fail unless there is both reference
        # and bipolar electrograms as well as ecgs. Annotations
        # are also required. And the gain of each signal.
        
        # Plot electrical data for this point
        current_index = self.annotate_dock.egm_selection.currentIndex()
        electric = self.system_manager.active_system.case.electric

        reference = electric.reference_egm.egm[current_index]
        bipolar = electric.bipolar_egm.egm[current_index]
        ecg = electric.ecg.ecg[current_index]
        signals = np.asarray([reference, bipolar, ecg])

        reference_gain = electric.reference_egm.gain[current_index]
        bipolar_gain = electric.bipolar_egm.gain[current_index]
        ecg_gain = electric.ecg.gain[current_index]
        signal_gains = np.asarray([reference_gain, bipolar_gain, ecg_gain])

        labels = np.asarray(["Ref", "Bipolar", "ECG"])
        for label, signal, gain in zip(labels, signals, signal_gains):
            artist = self.annotate_dock.signal_artists[label]
            artist._original_ydata = signal
            artist.set_ydata(artist._ystart + signal * np.exp(gain))

        # add the annotations
        annotations = electric.annotations
        times = self.egm_times

        local_annotation = annotations.local_activation_time[current_index]
        local_annotation_index = np.searchsorted(times, local_annotation)
        local_annotation_index = min(local_annotation_index, times.size)
        self.annotate_dock.update_annotation(
            signal=self.annotate_dock.signal_artists['Bipolar'],
            annotation=self.annotate_dock.annotation_artists['local_annotation_point'],
            annotation_line=self.annotate_dock.annotation_artists['local_annotation_line'],
            index=local_annotation_index,
        )

        reference_annotation = annotations.reference_activation_time[current_index]
        reference_annotation_index = np.searchsorted(times, reference_annotation)
        reference_annotation_index = min(reference_annotation_index, times.size)
        self.annotate_dock.update_annotation(
            signal=self.annotate_dock.signal_artists['Ref'],
            annotation=self.annotate_dock.annotation_artists['reference_annotation_point'],
            annotation_line=self.annotate_dock.annotation_artists['reference_annotation_line'],
            index=reference_annotation_index,
        )
        
        start_woi, stop_woi = annotations.window_of_interest[current_index]
        start_woi += reference_annotation
        stop_woi += reference_annotation
        self.annotate_dock.update_window_of_interest(start_woi, stop_woi)

        self.annotate_dock._update_scrollbar_from_woi(start_woi)
        self.annotate_dock.blit_artists()
        #self.annotate_dock.canvas.draw()
        #self.annotate_dock._initialise_scrollbar(start_woi=start_woi)
        #self.annotate_dock.update_active_artist()  # ensure the annotations/woi are added to the background

        # TODO: Display this point in the 3d viewer (e.g. render as large blue sphere)

    def annotation_on_button_press(self, event):
        """Update active artist on mouse button press, or set up call backs for moving annotations."""
        
        if event.button != MouseButton.LEFT:
            return
        
        annotater = self.annotate_dock
        
        if annotater.active_annotation_artist is not None:
            
            annotater.canvas.mpl_disconnect(annotater.cid_motion_notify_event_cursor_style)
            annotater.cid_motion_notify_event_annotation_position = annotater.canvas.mpl_connect(
                'motion_notify_event',
                self._annotation_on_mouse_move_position,
            )
            annotater.cid_button_release_event = annotater.canvas.mpl_connect(
                'button_release_event',
                self._annotation_on_button_release,
            )
            
            return
        
        # Don't do anything if the zoom/pan tools have been enabled.
        if annotater.canvas.widgetlock.locked():
            return
            
        if event.inaxes is None:
            return
        
        artist = annotater._get_signal_artist_under_point(
            cursor_position=np.asarray([[event.x, event.y]]),
        )
        
        # Don't do anything if the click is too far from an artist
        if artist is None:
            return

        annotater.active_signal_artist = artist.get_label()
        annotater.update_active_artist()

    def _annotation_on_mouse_move_position(self, event):
        """Update the active annotation line"""
        
        if event.inaxes is None:
            return
        
        artist_label = self.annotate_dock.active_annotation_artist
        artist = self.annotate_dock.annotation_artists[artist_label]
        new_time = int(event.xdata)
        
        if artist_label == "reference_annotation_line":
            
            current_time = np.asarray(artist.get_xdata())
            time_difference = new_time - current_time

            start_woi = self.annotate_dock.annotation_artists["start_woi"]
            start_woi.set_xdata(np.asarray(start_woi.get_xdata()) + time_difference)
            
            stop_woi = self.annotate_dock.annotation_artists["stop_woi"]
            stop_woi.set_xdata(np.asarray(stop_woi.get_xdata()) + time_difference)            

            self.annotate_dock.annotation_artists["reference_annotation_point"].set_xdata([new_time, new_time])
            self.annotate_dock._update_annotation_ydata(
                signal=self.annotate_dock.signal_artists["Ref"],
                annotation=self.annotate_dock.annotation_artists["reference_annotation_point"],
            )
        
        elif artist_label == "local_annotation_line":
            
            self.annotate_dock.annotation_artists["local_annotation_point"].set_xdata([new_time, new_time])
            self.annotate_dock._update_annotation_ydata(
                signal=self.annotate_dock.signal_artists["Bipolar"],
                annotation=self.annotate_dock.annotation_artists["local_annotation_point"],
            )

        artist.set_xdata([new_time, new_time])
        self.annotate_dock.blit_artists()

    def _annotation_on_button_release(self, event):
        """Disconnect callbacks for moving annotaitons, and set up callback for selecting active signal artist.
        
        Also store the annotation data in the active case.
        """
        
        annotater = self.annotate_dock
        annotater.canvas.mpl_disconnect(annotater.cid_motion_notify_event_annotation_position)
        annotater.cid_motion_notify_event_cursor_style = annotater.canvas.mpl_connect(
            'motion_notify_event',
            annotater._on_mouse_move_cursor_style,
        )
        annotater.canvas.mpl_disconnect(annotater.cid_button_release_event)
        
        # We also need to store the updated data
        current_index = self.annotate_dock.egm_selection.currentIndex()
        annotations = self.system_manager.active_system.case.electric.annotations
        annotation_artists = annotater.annotation_artists
        
        annotations.local_activation_time[current_index] = annotation_artists['local_annotation_point'].get_xdata()[0]
        ref_annotation = annotation_artists['reference_annotation_point'].get_xdata()[0]
        annotations.reference_activation_time[current_index] = ref_annotation
        woi = np.asarray([
            annotation_artists['start_woi'].get_xdata()[0] - ref_annotation,
            annotation_artists['stop_woi'].get_xdata()[0] - ref_annotation,
        ])
        annotations.window_of_interest[current_index] = woi
        
        # And re-interpolate onto the surface
        self.update_scalar_fields(event=None)

    def annotation_on_key_press(self, event):
        """Bindings for key press events in the annotation viewer"""
        
        sys.stdout.flush()
        
        if event.key not in [ 'down', 'up', 'r', 'w', 'W', 'l']:
            return

        # First check whether we need to change the electrogram that is displayed
        current_index = self.annotate_dock.egm_selection.currentIndex()

        if event.key == 'down':
            # move the the next electrogram
            n_items = self.annotate_dock.egm_selection.count()
            new_index = min(n_items-1, current_index+1)
            self.annotate_dock.egm_selection.setCurrentIndex(new_index)
            
        
        elif event.key == 'up':
            # move to the previous electrogram
            new_index = max(0, current_index - 1)
            self.annotate_dock.egm_selection.setCurrentIndex(new_index)

        # If we're chaning a point/line, ensure the key press was inside the figure
        if event.xdata is None:
            return
        
        time_index = np.searchsorted(self.egm_times, event.xdata)
        if time_index == self.egm_times.size:
            time_index -= 1
        time = self.egm_times[time_index]
        
        electric = self.system_manager.active_system.case.electric

        # TODO: Annotation times are signal indices. Collected at a specific frequency.
        #       Should be plotted a ms on the x-axis, taking into accoutnt the paper size
        #       and paper speed.
        if event.key == 'r':
            
            electric.annotations.reference_activation_time[current_index] = time_index
            
            voltage = electric.reference_egm.egm[current_index, time_index]  # + 2  # y offset
            gain = electric.reference_egm.gain[current_index]
            self.annotate_dock.update_reference_annotation(time, voltage, gain)
            self.annotate_dock.blit_artists()
            return

        elif event.key == 'w':
            
            # move the lower window boundary
            reference_annotation = electric.annotations.reference_activation_time[current_index]
            woi = electric.annotations.window_of_interest[current_index]
            woi[0] = time_index - reference_annotation
            woi = np.sort(woi)
            electric.annotations.window_of_interest[current_index] = woi
            self.annotate_dock.update_window_of_interest(*woi + reference_annotation)
            self.annotate_dock.blit_artists()
            return
        
        elif event.key == 'W':
            
            # move the upper boundary
            reference_annotation = electric.annotations.reference_activation_time[current_index]
            woi = electric.annotations.window_of_interest[current_index]
            woi[1] = time_index - reference_annotation
            woi = np.sort(woi)
            electric.annotations.window_of_interest[current_index] = woi
            self.annotate_dock.update_window_of_interest(*woi + reference_annotation)
            self.annotate_dock.blit_artists()
            return

        elif event.key == 'l':
            
            electric.annotations.local_activation_time[current_index] = time_index
            
            voltage = electric.bipolar_egm.egm[current_index, time_index]  # + 6  # y offset
            gain = electric.bipolar_egm.gain[current_index]
            self.annotate_dock.update_local_annotation(time, voltage, gain)
            self.annotate_dock.blit_artists()
            return

    def annotation_on_scroll_wheel(self, event):
        """Set the gain of the active line in the annotation viewer"""
        
        current_index = self.annotate_dock.egm_selection.currentIndex()
        electric = self.system_manager.active_system.case.electric
        gain_diff = 0.1 * (event.step * -1)  # scrolling up will decrease gain, down will increase gain
        
        if self.annotate_dock.active_signal_artist == "Ref":
            electric.reference_egm.gain[current_index] += gain_diff
            gain = electric.reference_egm.gain[current_index]
            
        elif self.annotate_dock.active_signal_artist == "Bipolar":
            electric.bipolar_egm.gain[current_index] += gain_diff
            gain = electric.bipolar_egm.gain[current_index]
            
        elif self.annotate_dock.active_signal_artist == "ECG":
            electric.ecg.gain[current_index] += gain_diff
            gain = electric.ecg.gain[current_index]
        
        self.annotate_dock.update_gain(gain)
        self.annotate_dock.blit_artists()

    def update_annotation_woi(self):
        """Set the woi based on the current value of the range slider"""

        # TODO: create a button/toolbar/menu option to update the woi/
        # reference annotation/LAT based on the current values

        current_index = self.annotate_dock.egm_selection.currentIndex()
        start_woi, stop_woi = self.annotate_dock.woi_slider.val

        annotations = self.system_manager.active_system.case.electric.annotations
        reference_annotation = annotations.reference_activation_time[current_index]
        
        annotations.window_of_interest[current_index, 0] = start_woi - reference_annotation
        annotations.window_of_interest[current_index, 1] = stop_woi - reference_annotation

    def _hide_mapping_points(self):
        """Show presiously deleted mapping points in the recycle bin, not the mapping points list"""

        non_recycled_rows = np.argwhere(self.mapping_points.model._include).squeeze()
        for row_number in non_recycled_rows:
            self.recycle_bin.table.setRowHidden(row_number, True)
        
        recycled_rows = np.argwhere(self.mapping_points.model._include == False).squeeze()
        for row_number in recycled_rows:
            self.mapping_points.table.setRowHidden(row_number, True)

    def delete_mapping_points(self, table, restore):
        """Move the selected mapping points into the recycle bin."""

        rows = table.selectionModel().selectedRows()
        row_numbers = [row.row() for row in rows]
        point_indices = np.asarray([row.sibling(row_number, 0).data() for row, row_number in zip(rows, row_numbers)], dtype=int)
        self.mapping_points.model._include[point_indices] = restore
        
        for row_number in row_numbers:
            self.mapping_points.table.setRowHidden(row_number, not restore)
            self.recycle_bin.table.setRowHidden(row_number, restore)

    def sort_mapping_points_and_recycle_bin(self, table, table_to_sort):
        """When one table is sorted by a column, sort the other table by the same column"""

        # prevent infinite recursion
        if table_to_sort._ignore_sort_signal:
            return
        
        table._ignore_sort_signal = True
        table_to_sort._ignore_sort_signal = True

        # sort second table based on first table
        sort_column = table.table.horizontalHeader().sortIndicatorSection()
        sort_order = table.table.horizontalHeader().sortIndicatorOrder()
        table_to_sort.table.sortByColumn(sort_column, sort_order)


        table._ignore_sort_signal = False
        table_to_sort._ignore_sort_signal = False

    def highlight_menubars_of_active_plotters(self):
        """Make the menubars of all docks in the active system blue."""

        for system_name, system in self.system_manager.systems.items():
            if system_name == self.system_manager.active_system.name:
                for dock in system.docks:
                    dock.main.menubar.setStyleSheet(
                        "QMenuBar{"
                        "background-color: #60798B;"
                        "}"
                    )
            else:
                for dock in system.docks:
                    dock.main.menubar.setStyleSheet(
                        "QMenuBar{"
                        "background-color: #455364;"  # default heading colour for qdarkstyle
                        "}"
                    )
            

def main():

    # Create an instance of Qapplication
    app = QtWidgets.QApplication(sys.argv)
    app.setWindowIcon(QtGui.QIcon(openep.view.static.LOGO))
    
    # This is necessary for vtk widgets to work with qdockwidgets on macOS
    # See: https://gitlab.kitware.com/vtk/vtk/-/issues/18454
    app.setStyle('Fusion')

    # setup stylesheet
    app.setStyleSheet(qdarkstyle.load_stylesheet(qt_api='pyside2'))

    # Create an instance of GUI
    window = OpenEPGUI()
    window.showMaximized()

    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
