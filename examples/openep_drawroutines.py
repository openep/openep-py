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

import sys

sys.path.append("/home/jra21/work/source/repos/opep")

from openep import io as openep_io
from openep import case as openep_case
from openep import mesh_routines as openep_mesh
from openep import case_routines as case_routines
from openep import draw_routines as draw
import numpy as np
import matplotlib.pyplot as plt


filename = "/home/jra21/work/source/repos/opep/examples/data/new_dataset_2.mat"

ep_case = openep_io.load_case(filename)

# DrawVoltage Map
hsurf = draw.draw_map(
    ep_case,
    volt="bip",
    freeboundary_color="black",
    freeboundary_width=5,
    cmap="jet_r",
    minval=0,
    maxval=2,
    volt_below_color="brown",
    volt_above_color="magenta",
    nan_color="gray",
    plot=True,
)

# GetAnatomicalStructure and DrawFreeBoundaries on mesh
# hsurf1 = draw.get_anatomical_structures(ep_case,plot=True)
