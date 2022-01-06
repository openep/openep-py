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
import openep.view.plotters_ui
import openep.view.system_manager_ui
import openep.view.system_ui

import openep.view.plotters
import openep.view.system_manager

import openep.view.static

class OpenEPGUI(QtWidgets.QMainWindow):

    def __init__(self):

        super().__init__()

        self._init_systems()
        self._init_ui()
        self._create_system_manager_ui()
        self._create_egm_canvas_dock()
        self._create_analysis_canvas_dock()
        self._add_dock_widgets()
        self._disable_dock_widgets()

    def _init_systems(self):
        """Containers and variables for keeping track of the systems loaded into the GUI."""

        self.system_manager = openep.view.system_manager.SystemManager()

        self.egm_slider = None
        self.egm_times = None
        self.egm_selection_filter = re.compile(r'\d+')

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

        self.system_manager_ui = openep.view.system_manager_ui.SystemManagetDockWidget(
            title="Systems"
        )

        # We also need to connect up the actions
        self.system_manager_ui.main.load_openep_mat_action.triggered.connect(self._load_openep_mat)
        self.system_manager_ui.main.load_opencarp_action.triggered.connect(self._load_opencarp)

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
        self.addDockWidget(Qt.RightDockWidgetArea, self.system_manager_ui)
        self.tabifyDockWidget(self.analysis_dock, self.system_manager_ui)

        for dock in [self.egm_dock, self.analysis_dock, self.system_manager_ui]:
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

    def _load_openep_mat(self):
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
        # Instead, this is done in self._add_data_to_openCARP
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
        dialogue.setWindowTitle('Add an openCARP data file')
        dialogue.setDirectory(QtCore.QDir.currentPath())
        dialogue.setFileMode(QtWidgets.QFileDialog.ExistingFile)
        dialogue.setNameFilter("openCARP data file (*.dat)")

        if dialogue.exec_():

            filename = dialogue.selectedFiles()[0]

            if data_type=="unipolar_egm":
                unipolar = np.loadtxt(filename)
                system.case.add_unipolar_electrograms(
                    unipolar=unipolar,
                    add_bipolar=True,
                    add_annotations=True,
                )

            system._determine_available_fields()  # ensure we can access the loaded data
        
            if len(system.plotters) == 0:

                # We need to dynamically add an option for creating 3D viewers to the main menubar
                add_view_for_system_action = QtWidgets.QAction(str(system.basename), self)
                add_view_for_system_action.triggered.connect(lambda: self.add_view(system))
                self.system_manager_ui.main.add_view_menu.addAction(add_view_for_system_action)

                # And also create the first 3d viewer
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
            case=system.case,
            prefix=prefix,
        )

    def update_system_manager_table(self, system):
        """Update the System manager table when a new system is loaded."""

        # Not sure why * 10 works. But if not multiplied by at least 10, there is vertical overlap between rows
        row_number = len(self.system_manager.systems) * 10

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

        try :
            self.system_manager.update_system_name(system=system, new_name=new_name)

        except KeyError as e:
       
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
            for dock, mesh  in zip(system.docks, system.meshes):
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
        _, *plotter_widgets = openep.view.plotters_ui.create_plotter_layout(
            plotter=plotter,
            add_link_views_button=plotter_is_secondary_view,
        )
        plotter.lower_limit, plotter.upper_limit, plotter.opacity, plotter.link_view_with_primary = \
            plotter_widgets

        dock = openep.view.system_ui.create_system_dock(plotter=plotter)

        # Add a Field menu to the menubar. This is used for selecting the scalar field to project onto the surface.
        dock, plotter = openep.view.system_ui.add_field_menu(
            dock=dock,
            plotter=plotter,
            system_name=system.name,
            scalar_fields=system.scalar_fields,
        )

        # Add a Show/Hide menu to the menubar. This is used for showing/hiding the surface mesh, mapping points,
        # and surface-projected mapping points
        dock, plotter = openep.view.system_ui.add_show_menu(
            dock=dock,
            plotter=plotter,
        )

        mesh = system.create_mesh()
        add_mesh_kws, add_points_kws = system._create_default_kws()
        free_boundaries = openep.mesh.get_free_boundaries(mesh)

        plotter.lower_limit.returnPressed.connect(lambda: self.update_colourbar_limits(system, index=index))
        plotter.upper_limit.returnPressed.connect(lambda: self.update_colourbar_limits(system, index=index))
        plotter.opacity.valueChanged.connect(lambda: self.update_opacity(system, index=index))

        # Link views across plotters
        if plotter_is_secondary_view:
            system.plotters[0].link_views_across_plotters(plotter)
            plotter.link_view_with_primary.toggled.connect(lambda: self.link_views_across_plotters(system, index=index))

        system.docks.append(dock)
        system.plotters.append(plotter)
        system.meshes.append(mesh)
        system.add_mesh_kws.append(add_mesh_kws)
        system.free_boundaries.append(free_boundaries)

        # We can't put this into a for loop.
        # For some reason, when iterating over the items in plotter.scalar_field_actions,
        # all actions are passed the key in the iteration for the scalars argument
        # TODO: Don't distinguish between Clinical and non-clinical data. Use clinical by default. If
        #       values are interpolated from mapping points onto the surface, then overwrite the
        #       default clinical values. Add an option to reset the interpolated values to the original
        #       clinical ones. This can be done by simply reading the file.
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
        # new point_data array as the active scalars. We need to undo this unwanted behaviour.
        mesh.set_active_scalars(active_scalars_name)

        self.draw_points(
            points=system.case.electric.bipolar_egm.points - system.case._mesh_center,
            plotter=plotter,
            actor_name="Mapping points",
            add_points_kws=add_points_kws,
        )
        plotter.renderer._actors['Mapping points'].SetVisibility(False)

        self.draw_points(
            points=system.case.electric.surface.nearest_point - system.case._mesh_center,
            plotter=plotter,
            actor_name="Surface-projected mapping points",
            add_points_kws=add_points_kws,
        )
        plotter.renderer._actors['Surface-projected mapping points'].SetVisibility(False)

        plotter.show_actions['Surface'].triggered.connect(
            lambda: self.update_actor_visibility(plotter, actor_name="Surface")
        )
        plotter.show_actions['Mapping points'].triggered.connect(
            lambda: self.update_actor_visibility(plotter, actor_name="Mapping points")
        )
        plotter.show_actions['Surface-projected mapping points'].triggered.connect(
            lambda: self.update_actor_visibility(plotter, actor_name="Surface-projected mapping points")
        )

        # If this is the first 3d viewer for the first system loaded, we need to update the egms etc.
        if (system.name == self.system_manager.active_system.name) and (len(system.plotters) == 1):
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
            if actor_name=="Surface" or actor_name.startswith('free_boundary_'):
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

        # TODO: submit a PR to pyvist to add a 'keep_camera_position' argument to 'plotter.unlink_views'
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

        self.remove_woi_slider()
        self.check_available_electrograms()
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

    def draw_points(self, points, plotter, actor_name, add_points_kws):
        """Render mapping points as 3D spheres.

        Args:
            points (np.ndarray): 3D coordinates of points to add
            plotter (BackgroundPlotter): plotter to which the points will be added 
            add_points_kws (dict): keyword arguments to pass to `plotter.add_points`
        """

        plotter.add_points(
            points,
            name=actor_name,
            **add_points_kws,
        )

    def check_available_electrograms(self):
        """Update the electrom types that can be plotted."""

        has_reference, has_bipolar, has_unipolar = self.system_manager.check_available_electrograms()

        self.has_electrograms = any([has_reference, has_bipolar, has_unipolar])
        self.egm_dock.setEnabled(self.has_electrograms)
        self.egm_reference_checkbox.setEnabled(has_reference)
        self.egm_bipolar_checkbox.setEnabled(has_bipolar)
        self.egm_unipolar_A_checkbox.setEnabled(has_unipolar)
        self.egm_unipolar_B_checkbox.setEnabled(has_unipolar)

        self._set_electrogram_times()

    def _set_electrogram_times(self):
        """Generate the electrogram times based on the """

        self.egm_times = self.system_manager.electrogram_times()

    def update_electrograms_from_text(self):
        """
        Extract electrograms at specified indices and re-plot.

        The points are specified by the user and set via use of the
        self.egm_select (QLineEdit) and button_egm_select (QPushButton) widgets.
        """

        # Get data for new set of points
        egm_text = self.egm_selection.text()
        selected_indices = self.egm_selection_filter.findall(egm_text)
        self.system_manager.active_system.egm_points = np.asarray(selected_indices, dtype=int)

        self.extract_electrograms()
        self.plot_electrograms()

    def update_electrograms_from_stored(self):

        # Retrieve data for the stored set of points
        egm_points = self.system_manager.active_system.egm_points
        egm_text = " ".join(egm_points.astype(str))
        self.egm_selection.setText(egm_text)

        self.extract_electrograms()
        self.plot_electrograms()

    def extract_electrograms(self):

        system = self.system_manager.active_system

        if self.egm_reference_checkbox.isEnabled():
            self.egm_reference_traces, self.egm_names = openep.case.get_electrograms_at_points(
                system.case,
                within_woi=False,
                buffer=0,
                indices=system.egm_points,
                egm_type="reference",
                return_names=True,
                return_lat=False,
            )

        if self.egm_bipolar_checkbox.isEnabled():
            self.egm_bipolar_traces, self.egm_names = openep.case.get_electrograms_at_points(
                system.case,
                within_woi=False,
                buffer=0,
                indices=system.egm_points,
                egm_type="bipolar",
                return_names=True,
                return_lat=False,
            )

        if self.egm_unipolar_A_checkbox.isEnabled():
            unipolar_traces, self.egm_names = openep.case.get_electrograms_at_points(
                system.case,
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
        separations = np.arange(self.system_manager.active_system.egm_points.size) * 2

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

        annotations = self.system_manager.active_system.case.electric.annotations

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

        annotations = self.system_manager.active_system.case.electric.annotations

        # from https://matplotlib.org/devdocs/gallery/widgets/range_slider.html
        start_woi, stop_woi = val
        self.egm_slider_lower_limit.set_xdata([start_woi, start_woi])
        self.egm_slider_upper_limit.set_xdata([stop_woi, stop_woi])

        reference_annotation = annotations.reference_activation_time[0]
        annotations.window_of_interest[:, 0] = start_woi - reference_annotation
        annotations.window_of_interest[:, 1] = stop_woi - reference_annotation

    def highlight_menubars_of_active_plotters(self):
        """Make the menubars of all docks in the active system blue."""

        for system_name, system in self.system_manager.systems.items():
            if system_name == self.system_manager.active_system.name:
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
    app.setWindowIcon(QtGui.QIcon(openep.view.static.LOGO))

    # Create an instance of GUI
    window = OpenEPGUI()
    window.showMaximized()

    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
