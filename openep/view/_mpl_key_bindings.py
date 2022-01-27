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
Disable the default matplotlib bindings to key press events.

If this is run as a script (`python _mpl_key_bindings`) it
will print all matplotlib bindings after removal.
"""

import matplotlib as mpl

def disable_all_bindings():

    mpl.rcParams['keymap.back'].remove('left')
    mpl.rcParams['keymap.back'].remove('c')
    mpl.rcParams['keymap.back'].remove('backspace')
    mpl.rcParams['keymap.back'].remove('MouseButton.BACK')

    mpl.rcParams['keymap.copy'].remove('ctrl+c')
    mpl.rcParams['keymap.copy'].remove('cmd+c')

    mpl.rcParams['keymap.forward'].remove('right')
    mpl.rcParams['keymap.forward'].remove('v')
    mpl.rcParams['keymap.forward'].remove('MouseButton.FORWARD')

    mpl.rcParams['keymap.fullscreen'].remove('f')
    mpl.rcParams['keymap.fullscreen'].remove('ctrl+f')

    #mpl.rcParams['keymap.grid'].remove('g')
    #mpl.rcParams['keymap.grid_minor'].remove('G')

    mpl.rcParams['keymap.help'].remove('f1')

    mpl.rcParams['keymap.home'].remove('h')
    mpl.rcParams['keymap.home'].remove('r')
    mpl.rcParams['keymap.home'].remove('home')
    
    mpl.rcParams['keymap.pan'].remove('p')

    mpl.rcParams['keymap.quit'].remove('ctrl+w')
    mpl.rcParams['keymap.quit'].remove('cmd+w')
    mpl.rcParams['keymap.quit'].remove('q')

    mpl.rcParams['keymap.save'].remove('s')
    mpl.rcParams['keymap.save'].remove('ctrl+s')

    mpl.rcParams['keymap.xscale'].remove('k')
    mpl.rcParams['keymap.xscale'].remove('L')

    mpl.rcParams['keymap.yscale'].remove('l')

    mpl.rcParams['keymap.zoom'].remove('o')

if __name__ == '__main__':

    disable_all_bindings()

    for k,v in mpl.rcParams.items():
        if -1 != k.find("keymap"):
            print("rcParams[%s]=%s"%(k,v))
