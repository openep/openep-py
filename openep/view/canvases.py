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
Create and manipulate Matplotlib canvases.
"""

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.pyplot as plt

import openep.view.custom_widgets

def create_canvas():
    
    figure, axis = plt.subplots(ncols=1, nrows=1)
    figure.set_facecolor("white")
    
    # hide it until we have data to plot
    axis.axis('off')
    # and don't display xy coordinates in the toolbar when hovering over the axis
    axis.format_coord = lambda x, y: ""

    canvas = FigureCanvas(figure)
    
    return canvas, figure, axis

def add_toolbar(canvas, parent, keep_actions=None):

    toolbar = openep.view.custom_widgets.CustomNavigationToolbar(
        canvas,
        parent,
        keep_actions=keep_actions,
    )

    return toolbar
