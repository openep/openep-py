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

"""Module containing functions to load an OpenEP dataset."""

import os
import scipy.io

from .matlab import _load_mat_v73, _load_mat_below_v73
from ..case.case import Case

__all__ = ["load_case", "load_mat"]


def _check_mat_version_73(filename):
    """Check if a MATLAB file is of version 7.3"""

    byte_stream, file_opened = scipy.io.matlab.mio._open_file(filename, appendmat=False)
    major_version, minor_version = scipy.io.matlab.miobase.get_matfile_version(byte_stream)

    return major_version == 2


def load_mat(filename):
    """Load a MATLAB file."""

    if _check_mat_version_73(filename):
        data = _load_mat_v73(filename)
    else:
        data = _load_mat_below_v73(filename)

    # These are indices
    data['surface']['triRep']['Triangulation'] -= 1

    return data


def load_case(filename, name=None):
    """
    Load a Case object from the given v7.3 MATLAB file.

    This assumes a number of names for objects found in the file.
    """
    data = load_mat(filename)

    if name is None:
        name = os.path.basename(filename)

    points = data['surface']['triRep']['X']
    indices = data['surface']['triRep']['Triangulation']
    act, bip = data['surface']['act_bip'].T
    uni, imp, frc = data['surface']['uni_imp_frc'].T

    # TODO: make classes for fields, electric, surface, rf, and rf_index
    fields = dict(act=act, bip=bip, uni=uni, imp=imp, frc=frc)
    electric = data['electric']
    surface = data['surface']

    try:
        rf = data['rf']
    except KeyError:
        rf = {}

    try:
        notes = data['notes']
    except KeyError:
        notes = []

    return Case(name, points, indices, fields, electric, surface, rf, notes)
