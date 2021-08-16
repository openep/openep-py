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

from openep import io as openep_io
from openep import case as openep_case
from openep import mesh_routines as openep_mesh
from openep import case_routines as case_routines
from openep import draw_routines as draw

import numpy as np

filename = '/home/jra21/work/source/repos/opep/examples/data/new_dataset_2.mat'
distance_thresh = 10

# Load opneep case
ep_case = openep_io.load_case(filename)

# Anatomic descriptions (Mesh) - nodes and indices
pts = ep_case.nodes
indices = ep_case.indices


# Electric data
# Locations – Cartesian co-ordinates, projected on to the surface 
locations = case_routines.get_electrogram_coordinates(ep_case,'type','bip')

i_egm = ep_case.electric['egm'].T
print('i_egm\n',i_egm)
i_vp = case_routines.getMappingPointsWithinWoI(ep_case)
print('i_vp\n',i_vp)
# macthing the shape of ivp with dataW
i_vp_egm = np.repeat(i_vp, repeats=i_egm.shape[1], axis=1)
print('i_vp_egm',i_vp_egm)
# macthing the shape of ivp with coords
i_vp_locations = np.repeat(i_vp, repeats=locations.shape[1],axis=1)

# Replacing the values outside the window of interest with Nan values
i_egm[~i_vp_egm] = np.nan
locations[~i_vp_locations] = np.nan

# For each mapping point, n, find the voltage amplitude
max_volt = np.amax(a=i_egm,axis=1).reshape(len(i_egm),1)
min_volt = np.amin(a=i_egm,axis=1).reshape(len(i_egm),1)

amplitude_volt = np.subtract(max_volt,min_volt)

for indx in range(amplitude_volt.shape[1]):
    temp_data = amplitude_volt[:,indx]
    temp_coords = locations
    i_nan = np.isnan(temp_data)
    temp_data=temp_data[~i_nan]
    temp_coords=temp_coords[~i_nan]


    interp = case_routines.OpenEPDataInterpolator(method='rbf',distanceThreshold=distance_thresh,rbfConstant=1)
    vertex_voltage_data = interp.interpolate(x0=temp_coords,d0=temp_data,x1=pts)



vsurf = draw.DrawMap(ep_case,
                    volt = vertex_voltage_data,
                    freeboundary_color='black',
                    cmap='jet_r',
                    freeboundary_width=5,
                    minval=0,
                    maxval=2,
                    volt_below_color='brown', 
                    volt_above_color='magenta', 
                    nan_color='gray', 
                    plot=True,)