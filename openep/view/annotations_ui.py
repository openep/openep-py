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
Create a dock widget for the annotation viewer.
"""
from re import A
from PySide2 import QtCore, QtWidgets

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backend_tools import Cursors
import matplotlib.pyplot as plt
from cycler import cycler
import mplcyberpunk

import numpy as np

from .custom_widgets import CustomDockWidget, CustomNavigationToolbar
from ._mpl_key_bindings import disable_all_bindings
from openep.case.case_routines import calculate_distance


plt.style.use("cyberpunk")

# mpl.rcParams['font.size'] = 8
# neon colors: teal/cyan, pink, yellow, matrix green, red, violet
# axes.prop_cycle: cycler('color', ['08F7FE', 'FE53BB', 'F5D300', '00ff41', 'r', '9467bd', ])
disable_all_bindings()


class AnnotationWidget(CustomDockWidget):
    """A dockable widget for annotating electrograms."""

    def __init__(
        self,
        title: str,
        active_linewidth: float,
        non_active_linewidth: float,
        woi_linewidth: float,
        annotation_linewidth: float,
        annotation_markersize: float,
        gain_prefactor: float,
    ):

        self._ACTIVE_LINEWIDTH = active_linewidth
        self._NON_ACTIVE_LINEWIDTH = non_active_linewidth
        self._WOI_LINEWIDTH = woi_linewidth
        self._ANNOTATION_LINEWITDH = annotation_linewidth
        self._ANNOTATION_MARKERSIZE = annotation_markersize
        self._GAIN_PREFACTOR = gain_prefactor
        self._SCROLLBAR_STEP = 0.1

        super().__init__(title)
        self._init_main_window()
        self._initialise_annotations()
        self._current_index = None  # index of mapping point currently plotted in the viewer

        # TODO: Get the resolution of the screen
        #       We can use this to set the size of the figure and the font size
        #       As well as the minimum size of the window        
        #app = QtWidgets.QApplication.instance()
        #screen_resolution = app.desktop().screenGeometry()
        #width, height = screen_resolution.width(), screen_resolution.height()

    def _init_main_window(self):
        
        self.main  = QtWidgets.QMainWindow()
        
        # The central widget will hold a matplotlib canvas and toolbar.
        # The canvas widget will also contain a QLabel that displays which
        # mapping point is currently selected
        self.canvas, self.figure, self.axes = self._init_canvas()
        self.signal_artists = {}  # dictionary of artists (lines, points, etc.)
        self.annotation_artists = {}  # dictionary of artists (for woi, annotations, etc.)
        self.active_signal_artist = None
        self.active_annotation_artist = None
        
        self.egm_selected, egm_selected_layout = self._init_selected()
        self.scrollbar = QtWidgets.QScrollBar(QtCore.Qt.Horizontal)
        self.scrollbar.actionTriggered.connect(self._update_limits_from_scrollbar)
        self.scrollbar_step = self._SCROLLBAR_STEP
        self.toolbar = CustomNavigationToolbar(
            canvas_=self.canvas,
            parent_=self,
            keep_actions=['Home', 'Zoom', 'Pan', 'Save'],
        )

        # Setting nested layouts
        central_widget = self._init_central_widget(
            egm_selected_layout, self.canvas, self.scrollbar, self.toolbar,
        )
        self.main.setCentralWidget(central_widget)
        self.setWidget(self.main)
        
        
        # Needed so we can detect key press events
        # See https://github.com/matplotlib/matplotlib/issues/707/#issuecomment-4181799
        self.canvas.setParent(central_widget)
        self.canvas.setFocusPolicy(QtCore.Qt.ClickFocus)
        self.canvas.setFocus()
        
    def _init_canvas(self):
        """
        Create an interactive matploblib canvas.
        """

        figure = plt.figure()
        figure.set_visible('off') # hide the figure until we have data to plot
        
        # ensure the axes takes up the entire figure
        axes = figure.add_axes([0,0,1,1])
        
        # only display x coordinate in the toolbar when hovering over the axis
        axes.format_coord = lambda x, y: f"{x:.1f} ms"

        canvas = FigureCanvas(figure)
        self.cid_draw_event = canvas.mpl_connect('draw_event', self._on_draw)

        return canvas, figure, axes

    def _init_selected(self):
        """Create a layout with widget that displays the point currently plotted.
        """

        annotate_selected = QtWidgets.QLabel()
        annotate_selected.setMinimumWidth(220)
        size_policy = QtWidgets.QSizePolicy()
        size_policy.setVerticalStretch(0)
        annotate_selected.setSizePolicy(size_policy)

        annotate_selected_layout = QtWidgets.QHBoxLayout()
        annotate_selected_layout.addWidget(annotate_selected)
        annotate_selected_layout.addStretch()
        
        return annotate_selected, annotate_selected_layout

    def _init_central_widget(self, egm_selected_layout, canvas, scrollbar, toolbar):
        """Create a placeholder widget to hold the toolbar and canvas.
        """
        
        central_layout = QtWidgets.QVBoxLayout()
        central_layout.addLayout(egm_selected_layout)
        central_layout.addWidget(canvas)
        central_layout.addWidget(scrollbar)
        central_layout.addWidget(toolbar)
        central_widget = QtWidgets.QWidget()
        central_widget.setLayout(central_layout)

        return central_widget

    def _initialise_annotations(self):
        """Create the annotations"""

        self._initialise_window_of_interest()
        self._initialise_reference_annotation()
        self._initialise_local_annotation()

    def _initialise_window_of_interest(self, start_woi=0, stop_woi=1):
        """Set default values for the sliders and plot the axvlines"""

        woi_slider_lower_limit = self.axes.axvline(
            start_woi,
            color='white',
            linestyle='--',
            linewidth=self._WOI_LINEWIDTH,
            alpha=0.6,
            label='start_woi',
            zorder=1,
            picker=False,
        )
        self.add_artist(
            artist=woi_slider_lower_limit,
            label="start_woi",
            signal=False,
        )        
        
        woi_slider_upper_limit = self.axes.axvline(
            stop_woi,
            color='white',
            linestyle='--',
            linewidth=self._WOI_LINEWIDTH,
            alpha=0.6,
            label='stop_woi',
            zorder=1,
            picker=False,
        )
        self.add_artist(
            artist=woi_slider_upper_limit,
            label="stop_woi",
            signal=False,
        )

    def _initialise_reference_annotation(self, time=0.5, voltage=6):
        """Plot a point at the reference activation time"""

        reference_annotation, = self.axes.plot(
            time,
            voltage,
            color='white',
            linewidth=0.0,
            marker='o',
            markersize=self._ANNOTATION_MARKERSIZE,
            label="reference_annotation_point",
            zorder=3,
            picker=False,
        )
        self.add_artist(
            artist=reference_annotation,
            label="reference_annotation_point",
            signal=False,
        )
        
        reference_annotation_line = self.axes.axvline(
            time,
            color='#08F7FE',
            linestyle='--',
            linewidth=self._ANNOTATION_LINEWITDH,
            alpha=0.6,
            label="reference_annotation_line",
            zorder=1,
            picker=False,
        )
        self.add_artist(
            artist=reference_annotation_line,
            label="reference_annotation_line",
            signal=False,
        )

    def _initialise_local_annotation(self, time=0.5, voltage=6):
        """Plot a point at the local activation time"""
        
        local_annotation, = self.axes.plot(
            time,
            voltage,
            color='white',
            linewidth=0,
            marker='o',
            markersize=self._ANNOTATION_MARKERSIZE,
            label="local_annotation_point",
            zorder=3,
            picker=False,
        )
        self.add_artist(
            artist=local_annotation,
            label="local_annotation_point",
            signal=False,
        )

        local_annotation_line = self.axes.axvline(
            time,
            color='#FE53BB',
            linestyle='--',
            linewidth=self._ANNOTATION_LINEWITDH,
            alpha=0.6,
            label="local_annotation_line",
            zorder=1,
            picker=False,
        )
        self.add_artist(
            artist=local_annotation_line,
            label="local_annotation_line",
            signal=False,
        )
    
    def _initialise_scrollbar(self, start_woi=None):
        """Set up the range and initial values of the scrollbar.
        
        Args:
            start_woi (int): Sample index at which woi start is located.
                If provided, this will be used to set the initial scrollbar
                position. If None, the initial scrollbar position will not
                be changed.
        """
        
        self.scrollbar_limits = np.asarray(self.axes.get_xlim()).astype(int)
        self.scrollbar.setPageStep(self.scrollbar_step * 100)
        
        if start_woi is not None:
            self._update_scrollbar_from_woi(start_woi)

        self._update_limits_from_scrollbar()

    def _update_scrollbar_from_woi(self, start_woi):
        """Move the scrollbar to be cover the woi lower limit."""
        
        lower_limit = start_woi - 30  # show the axes just to the left of the woi start
        slider_position = (lower_limit - self.scrollbar_limits[0]) / np.diff(self.scrollbar_limits)
        slider_index = int(slider_position * ((1 + self.scrollbar_step) * 100))
        slider_index = max(slider_index, 0)
        slider_index = min(slider_index, 99)
        self.scrollbar.setSliderPosition(slider_index)

    def _update_limits_from_scrollbar(self, event=None):
        """Update the displayed view of the scrollbar"""
        
        slider_position = self.scrollbar.sliderPosition() / ((1 + self.scrollbar_step) * 100)
        lower_limit = self.scrollbar_limits[0] + slider_position * np.diff(self.scrollbar_limits)
        upper_limit = lower_limit + np.diff(self.scrollbar_limits) * self.scrollbar_step
        self.axes.set_xlim(lower_limit, upper_limit)
        self.canvas.draw_idle()

    def _on_draw(self, event=None):
        """Store the background and blit other artists."""
        if event is not None:
            if event.canvas != self.canvas:
                raise RuntimeError
        self.background = self.canvas.copy_from_bbox(self.figure.bbox)
        self._draw_animated()

    def _draw_animated(self):
        """Draw all of the animated artists."""
        for artist in self.signal_artists.values():
            self.figure.draw_artist(artist)
        for artist in self.annotation_artists.values():
            self.figure.draw_artist(artist)

    def blit_artists(self):
        """Update the screen with animated artists."""

        canvas = self.canvas
        figure = self.figure

        # paranoia in case we missed the draw event,
        if self.background is None:
            self._on_draw(event=None)
        else:
            # restore the background
            canvas.restore_region(self.background)
            # draw all of the animated artists
            self._draw_animated()
            # update the GUI state
            canvas.blit(figure.bbox)

        # let the GUI event loop process anything it has to do
        canvas.flush_events()

    def _get_signal_artist_under_point(self, cursor_position):
        """Find the nearest artist to a given point.
        
        If it is within a given cutoff distance, return that artist.
        Otherwise, return None.
        """
        
        min_distance = 8
        nearest_artist = None
        
        for artist in self.signal_artists.values():
            
            artist_position = np.vstack(artist.get_data()).T  # axes coordinates, must transform to pixel coordinates
            distance = np.min(calculate_distance(cursor_position, self.axes.transData.transform(artist_position)))
            
            if distance < min_distance:
                min_distance = distance
                nearest_artist = artist
        
        return nearest_artist

    def _on_mouse_move_cursor_style(self, event):
        """Change the cursor style if the mouse moves over an annotation line"""
        
        # Don't do anything if the zoom/pan tools have been enabled.
        if self.canvas.widgetlock.locked():
            return
        
        artist = self._get_annotation_artist_under_point(
            cursor_position=np.asarray([[event.x, event.y]]),
        )
        
        if artist is None:
            self.active_annotation_artist = None
            self.canvas.set_cursor(Cursors.POINTER)
            return

        self.active_annotation_artist = artist.get_label()
        self.canvas.set_cursor(Cursors.RESIZE_HORIZONTAL)
    
    def _get_annotation_artist_under_point(self, cursor_position):
        """Find the nearest annotation artist to a given point"""
        
        min_distance = 6
        nearest_artist = None
        artists = [artist for artist_label, artist in self.annotation_artists.items() if not '_point' in artist_label]  # ignore points
    
        y0, y1 = self.axes.get_ylim()
        artist_y_position = np.linspace(y0, y1, 101)  # the artists are axvlines, so we need to generate the ydata

        for artist in artists:
            
            artist_x_position = np.repeat(artist.get_xdata()[0], repeats=artist_y_position.size)
            artist_position = np.vstack([artist_x_position, artist_y_position]).T  # axes coordinates, must transform to pixel coordinates
            distance = np.min(calculate_distance(cursor_position, self.axes.transData.transform(artist_position)))
            
            if distance < min_distance:
                min_distance = distance
                nearest_artist = artist

        return nearest_artist

    def add_artist(self, artist, label, signal=True):
        """Add an artist to be managed.
        
        Args:
            artist (Artists): Artist to be added. Will be set to 'animated'.
            label (str): Unique (and hashable) associated with the artist.
            signal (bool): Whether the artist is a signal (e.g. ecg, bipolar egm, etc.)
                or an annotation (e.g. window of interest, reference annotation, etc.).
                Defaults to True.
        """
        artist.set_animated(True)
        
        if signal:
            self.signal_artists[label] = artist
            return
            
        self.annotation_artists[label] = artist

    def update_active_artist(self):
        """Change the colour of the active artist"""
    
        for artist_label, artist in self.signal_artists.items():
            lw = self._ACTIVE_LINEWIDTH if artist_label == self.active_signal_artist else self._NON_ACTIVE_LINEWIDTH
            artist.set_linewidth(lw)

        self.blit_artists()

    def update_window_of_interest(self, start_woi, stop_woi):
        """Move vertical lines designating the window of interest"""

        self.annotation_artists['start_woi'].set_xdata([start_woi, start_woi])
        self.annotation_artists['stop_woi'].set_xdata([stop_woi, stop_woi])

    def update_reference_annotation(self, time, voltage, gain):
        """Plot the reference activation time"""
        
        ystart = self.signal_artists['Ref']._ystart
        scaled_voltage = ystart + np.exp(gain) * voltage
        self.annotation_artists['reference_annotation_point'].set_data([time], [scaled_voltage])
        self.annotation_artists['reference_annotation_line'].set_xdata([time, time])

    def update_local_annotation(self, time, voltage, gain):
        """Plot the local activation time"""
        
        ystart = self.signal_artists['Bipolar']._ystart
        scaled_voltage = ystart + np.exp(gain) * voltage
        self.annotation_artists['local_annotation_point'].set_data([time], [scaled_voltage])
        self.annotation_artists['local_annotation_line'].set_xdata([time, time])
    
    def update_annotation(self, signal, annotation, annotation_line, index):
        """Set the location of an annotation"""
        
        time = signal.get_xdata()[index]
        voltage = signal.get_ydata()[index]
        annotation.set_data([time], [voltage])
        annotation_line.set_xdata([time, time])

    def _update_annotation_ydata(self, signal, annotation):
        """After changing the gain of the signal, the ydata need to be modified"""
            
        timeseries, voltages = signal.get_data()
        activation_time = annotation.get_xdata()
        
        time_index = np.searchsorted(timeseries, activation_time)
        voltage = voltages[time_index]
        annotation.set_ydata(voltage)

    def update_gain(self, gain):
        """Set the gain of the active line"""

        # Update the signal
        label = self.active_signal_artist
        artist = self.signal_artists[label]
        original_ydata = artist._original_ydata
        artist.set_ydata(artist._ystart + original_ydata * np.exp(gain))
        
        # Update the y position of the annotation if necessary
        if label == "Ref":
            annotation_artist = self.annotation_artists['reference_annotation_point']
            self._update_annotation_ydata(signal=artist, annotation=annotation_artist)
        elif label == "Bipolar":
            annotation_artist = self.annotation_artists['local_annotation_point']
            self._update_annotation_ydata(signal=artist, annotation=annotation_artist)

    def plot_signals(self, times, signals, labels, gains):
        """Plot electrogram/ecg signals.
        
        Args:
            times (np.ndarray): Time at each point in the signal
            signals (np.npdarray): Signals to plot. 2D array of shape (n_signals, n_time_points)
            labels (np.ndarray): Name of each signal. We be used to label the y-axis.
        """
        
        self.axes.cla()
        self.axes.set_yticklabels([])
        self.axes.set_ylim(0, 12)

        # we need to horizontally shift the signals so they don't overlap
        y_start = 2
        y_separation = 4
        separations = y_start + np.arange(signals.shape[0]) * y_separation
        
        for signal, label, separation, gain in zip(signals, labels, separations, gains):
            
            line, = self.axes.plot(
                times,
                separation + np.exp(gain) * signal,
                linewidth=self._NON_ACTIVE_LINEWIDTH,
                label=label,
                zorder=2,
                picker=False,
                alpha=1,
            )
            
            # store the artists so we can blit later on
            line._original_ydata = signal.copy()
            line._ystart = separation
            self.add_artist(artist=line, label=label, signal=True)

        # Set an active artists (has a thicker linewidth than other lines. The gain can be set by scrolling).
        self.active_signal_artist = labels[0]

        self.axes.set_yticks(separations, labels)  # need to explicitly set the grid positions
        legend = self.figure.legend(
            handles=[line for line in self.signal_artists.values()],
            loc='upper center',
            ncol=3,
            frameon=False,
        )
        self._legend_handles = {handle.get_label(): handle for handle in legend.legendHandles}
        for handle in self._legend_handles.values():
            handle.set_linewidth(2.5)

        # Remove the border and ticks
        plt.tick_params(axis='both', which='both', length=0)
        for spine in ['left', 'right', 'top']:
            self.axes.spines[spine].set_visible(False)
        self.axes.spines['bottom'].set_alpha(0.4)
        self.axes.tick_params(axis=u'both', which=u'both',length=0)
        self.axes.grid(axis='y', visible=True, linestyle="-")
        self.axes.grid(axis='x', visible=False)
        self.axes.margins(0)  # don't add padding the the x axis when plotting

    def activate_figure(self, xmin, xmax):
        """Show the figure"""
        
        self.figure.set_visible(True)
        self.axes.set_xlim(xmin, xmax)
        self.axes.set_ylim(0, 12)
        self.canvas.draw()  # ensure the figure background is stored in self.background
    
    def deactivate_figure(self):
        """Clear the axes and hide the figure"""
        
        self.figure.set_visible(False)
