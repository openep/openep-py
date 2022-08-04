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

carp = openep.load_opencarp(
    points="/Users/paul/github/openep-misc/examples/data/pig21_endo_coarse.pts",
    indices="/Users/paul/github/openep-misc/examples/data/pig21_endo_coarse.elem",
    scale=1000,
)
unipolar = np.loadtxt("/Users/paul/github/openep-misc/examples/data/pig21_endo_coarse_phie.dat")
carp.add_unipolar_electrograms(unipolar=unipolar)

openep.export_openep_mat(carp, 'test-export.mat')
new_carp = openep.load_openep_mat('test-export.mat')