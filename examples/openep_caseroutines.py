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

import numpy as np
import openep


filename = "/Users/paul/github/openep-py/examples/data/new_dataset_2.mat"
case = openep.load_case(filename)
mesh = case.create_mesh()

# Get mapping points within WOI
mapping_points = openep.case.get_mapping_points_within_woi(case, indices=np.arange(10))

# Get electrograms
electrograms, names, local_activation_times = openep.case.get_electrograms_at_points(case, indices=1)

# Get interpolated voltages
interpolated_voltages = openep.case.get_voltage_electroanatomic(case)
