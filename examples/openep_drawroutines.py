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


import openep


filename = "/Users/paul/github/openep-py/examples/data/new_dataset_2.mat"
case = openep.load_case(filename)
mesh = case.create_mesh()

# DrawVoltage Map
plotter = openep.draw.draw_map(
    mesh=mesh,
    field=case.fields['bip'],
)
plotter.show()
