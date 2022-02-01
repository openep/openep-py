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
from PySide2 import QtCore, QtGui, QtWidgets
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backend_bases import MouseButton
import numpy as np

from .custom_widgets import CustomDockWidget, CustomNavigationToolbar
from ._mpl_key_bindings import disable_all_bindings


mpl.rcParams['font.size'] = 9
disable_all_bindings()

class AnnotationWidget(CustomDockWidget):
    """A dockable widget for annotating electrograms."""

    def __init__(self, title: str):

        super().__init__(title)
        self._init_main_window()
        self._initialise_window_of_interest()
        self._initialise_reference_annotation()
        self._initialise_local_annotation()
    
    def _init_main_window(self):
        
        self.main  = QtWidgets.QMainWindow()
        
        # The dock is set to have bold font (so the title stands out)
        # But all other widgets should have normal weighted font
        main_font = QtGui.QFont()
        main_font.setBold(False)
        self.main.setFont(main_font)
        
        # The central widget will hold a matplotlib canvas and toolbar.
        # The canvas widget will also contain a QComboBox for selecting
        # the electrogram to annnotate.
        self.canvas, self.figure, self.axes = self._init_canvas()
        self.artists = {}  # dictionary of artists (lines, points, etc.)
        self.active_artist_label = None
        
        self.egm_selection, egm_selection_layout = self._init_selection()
        canvas_layout = QtWidgets.QVBoxLayout(self.canvas)
        canvas_layout.addLayout(egm_selection_layout)
        canvas_layout.addStretch()

        toolbar = CustomNavigationToolbar(
            canvas_=self.canvas,
            parent_=self,
            keep_actions=['Home', 'Zoom', 'Pan', 'Save'],
        )

        # Setting nested layouts
        central_widget = self._init_central_widget(self.canvas, toolbar)
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

        figure, axes = plt.subplots(ncols=1, nrows=1)
        #figure.set_facecolor("white")
        figure.set_visible('off') # hide the figure until we have data to plot
        
        # only display x coordinate in the toolbar when hovering over the axis
        axes.format_coord = lambda x, y: f"{x:.1f} ms"

        canvas = FigureCanvas(figure)
        canvas.mpl_connect('draw_event', self._on_draw)
        canvas.mpl_connect('pick_event', self._on_pick)

        return canvas, figure, axes

    def _init_selection(self):
        """Create a layout with widgets for selecting which electrogram to annotate.
        """

        annotate_selection = QtWidgets.QComboBox()
        annotate_selection.setMinimumWidth(220)
        #annotate_selection.setStyleSheet('selection-background-color: red')
        #annotate_selection.setStyleSheet('border: 1px solid #d8dcd6; background-color: white;')

        annotate_selection.setStyleSheet(
            "QWidget{"
            "background-color: white;"
            "selection-background-color: #168CFF;"
            "border: 1px solid #d8dcd6;"
            "}"
        )


        annotate_selection_layout = QtWidgets.QHBoxLayout()
        annotate_selection_layout.addWidget(annotate_selection)
        annotate_selection_layout.addStretch()

        return annotate_selection, annotate_selection_layout

    def _init_central_widget(self, canvas, toolbar):
        """Create a placeholder widget to hold the toolbar and canvas.
        """
        
        central_layout = QtWidgets.QVBoxLayout()
        central_layout.addWidget(canvas)
        central_layout.addWidget(toolbar)
        central_widget = QtWidgets.QWidget()
        central_widget.setLayout(central_layout)
        central_widget.setStyleSheet("border-width: 0px; border: 0px; background-color: #d8dcd6;")
        
        return central_widget
    
    def _on_draw(self, event):
        """Store the background and blit other artists."""
        
        self.background = self.canvas.copy_from_bbox(self.axes.bbox)
        self._blit_artists()

    def _blit_artists(self):
        """Reload the stored background and blit the artists"""
        
        self.canvas.restore_region(self.background)
        for artist in self.artists.values():            
            self.axes.draw_artist(artist)
        
        self.axes.draw_artist(self.reference_annotation)
        self.axes.draw_artist(self._reference_annotation_line)
        self.axes.draw_artist(self.local_annotation)
        self.axes.draw_artist(self._local_annotation_line)
        
        self.canvas.blit(self.axes.bbox)
    
    def update_gain(self, gain):
        """Set the gain of the active line"""

        label = self.active_artist_label
        artist = self.artists[label]
        original_ydata = artist._original_ydata
        
        artist.set_ydata(artist._ystart + original_ydata * np.exp(gain))
        
        if label == "Ref":
            self._update_annotation_ydata(signal=artist, annotation=self.reference_annotation)
        elif label == "Bipolar":
            self._update_annotation_ydata(signal=artist, annotation=self.local_annotation)
        
        self._blit_artists()
    
    def _on_pick(self, event):
        """Set the active artist"""
        
        label = event.artist.get_label()
        if self.active_artist_label == label:
            return
        self.active_artist_label = label
        
        
        if event.artist.get_label() in ['woi_start', 'woi_stop']:
            self._update_window_of_interest(event)
            return
        
        self._update_active_artist()
    
    def _initialise_window_of_interest(self, start_woi=0, stop_woi=1):
        """Set default values for the sliders and plot the axvlines"""

        woi_slider_lower_limit = self.axes.axvline(
            start_woi,
            color='grey',
            linestyle='--',
            linewidth=0.6,
            alpha=0.6,
            label='woi_start',
            zorder=1,
            picker=True,
        )
        woi_slider_lower_limit.set_animated(True)
        self.artists['woi_start'] = woi_slider_lower_limit
        
        
        woi_slider_upper_limit = self.axes.axvline(
            stop_woi,
            color='grey',
            linestyle='--',
            linewidth=0.6,
            alpha=0.6,
            label='woi_stop',
            zorder=1,
            picker=True,
        )
        woi_slider_upper_limit.set_animated(True)
        self.artists['woi_stop'] = woi_slider_upper_limit
    
    def update_window_of_interest(self, start_woi, stop_woi):
        """Plot vertical lines designating the window of interest"""

        self.artists['woi_start'].set_xdata([start_woi, start_woi])
        self.artists['woi_stop'].set_xdata([stop_woi, stop_woi])
        self._blit_artists()

    def _update_window_of_interest(self, event):
        """This is called when the woi line is picked"""
        pass

    def _update_active_artist(self):
        """Change the colour of the active artist"""
        
        for artist_label, artist in self.artists.items():
            colour = 'xkcd:azure' if artist_label == self.active_artist_label else 'xkcd:steel blue'
            artist.set_color(colour)
        
        self._blit_artists()
        
    def _initialise_reference_annotation(self, time=0.5, voltage=6):
        """Plot a point at the reference activation time"""
        
        self.reference_annotation, = self.axes.plot(
            time,
            voltage,
            color='red',
            linewidth=0,
            marker='o',
            markersize=4,
            zorder=3,
            picker=False,
        )
        self.reference_annotation.set_animated(True)
        
        self._reference_annotation_line = self.axes.axvline(
            time,
            color='red',
            linestyle='--',
            linewidth=0.6,
            alpha=0.6,
            zorder=1,
            picker=False,
        )
        self._reference_annotation_line.set_animated(True)
    
    def update_reference_annotation(self, time, voltage, gain):
        """Plot the reference activation time"""
        
        ystart = self.artists['Ref']._ystart
        scaled_voltage = ystart + np.exp(gain) * voltage
        self.reference_annotation.set_data([time], [scaled_voltage])
        self._reference_annotation_line.set_xdata([time, time])
        self._blit_artists()
        
    def _initialise_local_annotation(self, time=0.5, voltage=6):
        """Plot a point at the local activation time"""
        
        self.local_annotation, = self.axes.plot(
            time,
            voltage,
            color='green',
            linewidth=0,
            marker='o',
            markersize=4,
            zorder=3,
            picker=False,
        )
        self.local_annotation.set_animated(True)

        self._local_annotation_line = self.axes.axvline(
            time,
            color='green',
            linestyle='--',
            linewidth=0.6,
            alpha=0.6,
            zorder=1,
            picker=False,
        )
        self._local_annotation_line.set_animated(True)

    def update_local_annotation(self, time, voltage, gain):
        """Plot the local activation time"""
        
        ystart = self.artists['Bipolar']._ystart
        scaled_voltage = ystart + np.exp(gain) * voltage
        self.local_annotation.set_data([time], [scaled_voltage])
        self._local_annotation_line.set_xdata([time, time])
        self._blit_artists()
    
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

    def initialise_egm_selection(self, selections):
        """Set the selections available in the QComboBox"""
        
        self.egm_selection.clear()
        self.egm_selection.insertItems(0, selections)

    def plot_signals(self, times, signals, labels, gains):
        """Plot electrogram/ecg signals.
        
        Args:
            times (np.ndarray): Time at each point in the signal
            signals (np.npdarray): Signals to plot. 2D array of shape (n_signals, n_time_points)
            labels (np.ndarray): Name of each signal. We be used to label the y-axis.
        """
        
        # First we need to clear the axis, but we want to keep to previous limits
        xlim = self.axes.get_xlim()
        ylim = self.axes.get_ylim()
        
        self.axes.cla()
        self.axes.set_yticklabels([])
        self.axes.set_ylim(ylim)
        self.axes.set_xlim(xlim)
        
        # we need to horizontally shift the signals so they don't overlap
        y_start = 2
        y_separation = 4
        separations = y_start + np.arange(signals.shape[0]) * y_separation
        
        for signal, label, separation, gain in zip(signals, labels, separations, gains):
            
            line, = self.axes.plot(
                times,
                separation + np.exp(gain) * signal,
                linewidth=0.8,
                label=label,
                zorder=2,
                picker=True,
                alpha=1,
            )
            
            # store the artists so we can blit later on
            self.artists[label] = line
            line._original_ydata = signal.copy()
            line._ystart = separation
            line.set_animated(True)

        # Set an active artists (has a different colour to to others. The gain can be set by scrolling).
        self.active_artist_label = labels[0]
        self.artists[self.active_artist_label].set_color('xkcd:azure')

        self.axes.set_yticks(separations, labels)

        # Remove the border and ticks
        plt.tick_params(axis='both', which='both', length=0)
        for spine in ['left', 'right', 'top']:
            self.axes.spines[spine].set_visible(False)
        self.axes.spines['bottom'].set_alpha(0.4)
        self.axes.tick_params(axis=u'both', which=u'both',length=0)
        self.axes.grid(axis='y')

    def activate_figure(self, xmin, xmax):
        """Show the figure"""
        
        self.figure.set_visible(True)
        self.axes.set_xlim(xmin, xmax)
        self.axes.set_ylim(0, 12)
        self.canvas.draw()  # ensure the figure background is stored in self.background
    
    def deactivate_figure(self):
        """Clear the axes and hide the figure"""
        
        self.figure.set_visible(False)
