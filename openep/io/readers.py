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
Loading datasets --- :mod:`openep.io.readers`
=============================================

This module contains functions to load an OpenEP dataset.

Note
----

`openep-py` is currently only able to load:
    * data exported from the MATLAB implementation of OpenEP. Further, the
     `rfindex` data is not yet loaded.
    * data from an openCARP simulation


Example of loading a dataset
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Case data can be loaded as follows:

.. code:: python

    import openep
    from openep._datasets.openep_datasets import DATASET_2_V73

    case = openep.load_openep_mat(DATASET_2_V73)

This will load the dataset and store the information in a `Case` object.
See :class:`openep.data_structures.case.Case` for information on the attributes
and methods of `Case`.

.. autofunction:: load_openep_mat
.. autofunction:: load_opencarp

"""

import os
import scipy.io

import numpy as np

from .matlab import _load_mat_v73, _load_mat_below_v73
from ..data_structures.surface import extract_surface_data, empty_fields
from ..data_structures.electric import extract_electric_data, empty_electric
from ..data_structures.ablation import extract_ablation_data, empty_ablation
from ..data_structures.case import Case

__all__ = ["load_openep_mat", "load_mat", "load_opencarp"]


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


def load_openep_mat(filename, name=None):
    """
    Load a Case object from a MATLAB file.

    Currently, cases can only be loaded from files created using the MATLAB
    implementation of OpenEP.

    Args:
        filename (str): path to MATLAB file to be loaded (including the .mat
            extension.)
        name (str): name to give this dataset. The default is `None`, in which case
            the filename is used at the name.
    
    Returns:
        case (Case): an OpenEP Case object that contains the surface, electric and
            ablation data.
    """
    data = load_mat(filename)

    if name is None:
        name = os.path.basename(filename)

    points, indices, fields = extract_surface_data(data['surface'])
    electric = extract_electric_data(data['electric'])

    try:
        ablation = extract_ablation_data(data['rf'])
    except KeyError:
        ablation = None

    try:
        notes = data['notes']
    except KeyError:
        notes = []

    return Case(name, points, indices, fields, electric, ablation, notes)


def load_opencarp(
    points,
    indices,
    name=None,
):
    """
    Load data from an OpenCARP simulation.

    Args:
        points (str): Path to the openCARP points file.
        indices (str): Path to the openCARP element file. Currently, only triangular meshes are
            supported.
        name (str, optional): Name of the dataset. If None, the basename of the points file
            will be used as the name.

    Returns:
        case (Case): an OpenEP Case object that contains the points and indices.

    Note
    ----
    All other attributes of the Case object will be empty numpy arrays.

    """

    name = os.path.basename(points) if name is None else name

    points_data = np.loadtxt(points, skiprows=1)
    indices_data = np.loadtxt(indices, skiprows=1, usecols=[1, 2, 3], dtype=int)  # ignore the tag for now

    fields = empty_fields()
    electric = empty_electric()
    ablation = empty_ablation()
    notes = []

    return Case(name, points_data, indices_data, fields, electric, ablation, notes)
