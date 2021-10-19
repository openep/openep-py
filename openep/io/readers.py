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
import numpy as np
import scipy.io
import h5py

from ..case.case import Case

__all__ = ["load_case"]


def _mat_version_supported(filename):

    byte_stream, file_opened = scipy.io.matlab.mio._open_file(filename, appendmat=False)
    major_version, minor_version = scipy.io.matlab.miobase.get_matfile_version(byte_stream)

    return major_version == 2


def _dereference_strings(file_pointer, references):
    """Resolve an array of references that point to strings."""

    return np.asarray([file_pointer.get(ref)[:].tobytes().decode('utf-16') for ref in references])


def _decode_string(ints):
    """Decode an array of unsigned integers to a single string.
    """

    return ints.tobytes().decode('utf-16').strip('\x00')


def _dereference_nans(file_pointer, references):
    """Resolve an array of references that point to nans/infs."""

    return np.asarray([file_pointer.get(ref)[:] for ref in references]).ravel()


def load_mat(filename):
    """
    Load a v7.3 MATLAB file.

    h5py is used to read the file.

    Currently, all references in the HDF5 file are resolved except for 'userdata/rfindex/grid'
    """

    if not _mat_version_supported(filename):
        raise NotImplementedError("Only MATLAB v7.3 files are currently supported.")

    strings_to_dereference = {
        'userdata/notes',
        'userdata/electric/electrodeNames_bip',
        'userdata/electric/electrodeNames_uni',
        'userdata/electric/names',
        'userdata/electric/tags',
    }

    strings_to_decode = {
        'userdata/cartoFolder',
        'userdata/rfindex/tag/index/name',
    }

    nans_to_dereference = {
        'userdata/electric/impedances/time',
        'userdata/electric/impedances/value',
    }

    dat = {}
    with h5py.File(filename, "r") as f:

        def _visitor(key, value):
            if (
                isinstance(value, h5py.Dataset)
                and not key.startswith('#refs#')
            ):

                values = value[:]

                if key in strings_to_dereference:
                    values = _dereference_strings(
                        file_pointer=f,
                        references=values.ravel(),
                    )

                elif key in strings_to_decode:
                    values = _decode_string(values.ravel())

                elif key in nans_to_dereference:
                    values = _dereference_nans(
                        file_pointer=f,
                        references=values.ravel(),
                    )

                # ignore references for now
                dat[key] = values

        f.visititems(_visitor)  # visit all items in the file and populate dat

    # References in the grid are not currently resolved
    dat.pop('userdata/rfindex/grid', None)

    return dat


def load_case(filename, name=None):
    """
    Load a Case object from the given v7.3 MATLAB file.

    This assumes a number of names for objects found in the file.
    """
    dat = load_mat(filename)

    if name is None:
        name = os.path.basename(filename)

    nodes = dat["userdata/surface/triRep/X"].T
    inds = dat["userdata/surface/triRep/Triangulation"].T - 1
    act, bip = dat["userdata/surface/act_bip"]
    uni, imp, frc = dat["userdata/surface/uni_imp_frc"]

    fields = dict(act=act, bip=bip, uni=uni, imp=imp, frc=frc)
    electric = {}
    rf = {}
    other_data = {}

    for k, v in dat.items():
        _, section, *objname = k.split("/", 2)

        if objname and section in ("electric", "rf"):
            dobj = rf if section == "rf" else electric
            dobj[objname[0]] = v
        else:
            other_data[k] = v

    return Case(name, nodes, inds, fields, electric, rf, other_data)
