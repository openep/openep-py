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
Create and manipulate PyVista-Qt BackgroundPlotters.
"""

from PyQt5 import QtCore, QtWidgets
from pyvistaqt import BackgroundPlotter


def create_plotter():
    """Initialise a BackgroundPlotter"""

    plotter = BackgroundPlotter(
        show=False,
        app=QtWidgets.QApplication.instance(),
        allow_quit_keypress=False,
        line_smoothing=False,
        point_smoothing=False,
        polygon_smoothing=False,
        toolbar=False,
        editor=False,
        menu_bar=False,
        title="Voltage",
        border=False,
        update_app_icon=False,
    )
    plotter.background_color = 'white'
    plotter.setMinimumSize(QtCore.QSize(50, 50))

    return plotter

def create_colourbar_widgets(plotter, lower, upper):
    """Create widgets for setting the colourbar limits of a BackgroundPlotter"""

    lower_limit = QtWidgets.QLineEdit(plotter)
    lower_limit.setStyleSheet("background-color: white; border: 1px solid lightGray;")
    lower_limit.setGeometry(200, 10, 50, 40)
    lower_limit.setText(str(lower))

    upper_limit = QtWidgets.QLineEdit(plotter)
    upper_limit.setStyleSheet("background-color: white; border: 1px solid lightGray;")
    upper_limit.setGeometry(260, 10, 50, 40)
    upper_limit.setText(str(upper))

    button_set_thresholds = QtWidgets.QPushButton("Set colourbar limits:", plotter)
    button_set_thresholds.setStyleSheet("background-color: lightGray")
    button_set_thresholds.setGeometry(10, 10, 180, 40)

    return lower_limit, upper_limit, button_set_thresholds
