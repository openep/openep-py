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
import matplotlib.pyplot as plt
import openep


filename = "/Users/paul/github/openep-py/examples/data/new_dataset_1.mat"
case = openep.load_case(filename)
mesh = case.create_mesh()

# Get mapping points within WOI
mapping_points = openep.case.get_mapping_points_within_woi(case, indices=np.arange(10))

# Get electrograms
electrograms, names, local_activation_times = openep.case.get_electrograms_at_points(
    case,
    indices=[1, 10, 100],
    return_names=True,
    return_lat=True,
)

# Plot the electrogram traces
# Array of times that are within the window of interest:
times = openep.case.get_woi_times(case)
# Array of times that are within the window of interest. Pacing at time = 0 seconds.
relative_times = openep.case.get_woi_times(case, relative=True)
# Now plot the traces
fig, axis = openep.draw.plot_electrograms(relative_times, electrograms[:, times], names=names)
plt.show()

# Plot Carto bipolar voltages
plotter = openep.draw.draw_map(
    mesh=mesh,
    field=case.fields.bipolar_voltage,
)
plotter.background_color = "white"
plotter.show()

# Interpolate bipolar voltage onto surface from raw electrograms
interpolated_voltages = openep.case.interpolate_voltage_onto_surface(case, max_distance=15)
plotter = openep.draw.draw_map(
    mesh=mesh,
    field=interpolated_voltages,
)
plotter.background_color = "white"
plotter.show()
