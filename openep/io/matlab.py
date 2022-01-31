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

from collections import defaultdict

import h5py
import scipy.io
import numpy as np

__all__ = []


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


def _visit_mat_v73_(file_pointer):
    """Extract all arrays from a HDF5 matlab file."""

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

    data = {}

    def _visitor(key, value):

        if (
            isinstance(value, h5py.Dataset)
            and not key.startswith('#refs#')
        ):

            values = value[:]

            if key in strings_to_dereference:
                values = _dereference_strings(
                    file_pointer=file_pointer,
                    references=values.ravel(),
                )

            elif key in strings_to_decode:
                values = _decode_string(values.ravel())

            elif key in nans_to_dereference:
                values = _dereference_nans(
                    file_pointer=file_pointer,
                    references=values.ravel(),
                )

            data[key] = values

    file_pointer.visititems(_visitor)  # visit all items in the file and populate dat

    return data


def _mat_v73_transform_arrays(data):
    """Flatten or transpose arrays if necessary. Cast arrays from float to int if necessary."""

    cast_to_int = {
        'userdata/electric/annotations/mapAnnot',
        'userdata/electric/annotations/referenceAnnot',
        'userdata/electric/annotations/woi',
        'userdata/surface/isVertexAtRim',
        'userdata/surface/triRep/Triangulation',
        'userdata/surface/triRep/Connectivity',
    }

    for key in data:

        if not isinstance(data[key], np.ndarray):
            continue

        if data[key].shape[0] == 1:
            data[key] = data[key].ravel()
        else:
            data[key] = data[key].T

        if key in cast_to_int:
            data[key] = data[key].astype(int)

        if key == 'userdata/electric/electrodeNames_uni':
            data[key] = data[key].reshape(2, -1).T

    return data


def _mat_v73_flat_to_nested(data):
    """Make a flat dictionary into a nested one.

    The forward slashes in the keys of the flat dictionary will be used to define the
    nesting points for keys in the nested dictionary. e.g. ['userdata/electric'] would become
    ['userdata']['electric].

    Args:
        data ([type]): [description]
    """

    nested_dict = lambda: defaultdict(nested_dict)  # noqa: E731

    nested_data = nested_dict()
    for key in data:

        nested_keys = key.split('/')[1:]

        # can this block be replaced by a recursive function?
        if len(nested_keys) == 1:
            key1 = nested_keys[0]
            nested_data[key1] = data[key]

        elif len(nested_keys) == 2:
            key1, key2 = nested_keys
            nested_data[key1][key2] = data[key]

        elif len(nested_keys) == 3:
            key1, key2, key3 = nested_keys
            nested_data[key1][key2][key3] = data[key]

        elif len(nested_keys) == 4:
            key1, key2, key3, key4 = nested_keys
            nested_data[key1][key2][key3][key4] = data[key]

        else:
            raise ValueError(f"Need more leaves! {key}")

    return nested_data


def _load_mat_v73(filename):
    """
    Load a v7.3 MATLAB file.

    h5py is used to read the file.

    Currently, all references in the HDF5 file are resolved except for 'userdata/rfindex/grid'
    """

    with h5py.File(filename, "r") as f:
        data = _visit_mat_v73_(f)

    # rfindex is a matlab class - not readable with Python
    data.pop('userdata/rfindex/tag', None)
    data.pop('userdata/rfindex/grid', None)

    # Some arrays need to be flattened or transposed
    data = _mat_v73_transform_arrays(data)

    # make the flat dictionary nested
    data = _mat_v73_flat_to_nested(data)

    return data


def _decode_empty_strings_array(arr):
    """"""
    return np.asarray([item.tobytes().decode('utf-16') for item in arr])


def _cast_to_float(arr):
    """Cast a numpy array to float"""

    return arr.astype(float)


def _load_mat_below_v73(filename):
    """
    Load a MATLAB file of version less than v7.3

    scipy.io.loadmat is used to read the file.
    """

    data = scipy.io.loadmat(
        filename,
        appendmat=False,
        mat_dtype=False,
        chars_as_strings=True,
        struct_as_record=False,
        squeeze_me=True,
        simplify_cells=True,
    )['userdata']

    data['electric']['tags'] = _decode_empty_strings_array(data['electric']['tags'])

    data['electric']['impedances']['time'] = _cast_to_float(data['electric']['impedances']['time'])
    data['electric']['impedances']['value'] = _cast_to_float(data['electric']['impedances']['value'])

    try:
        # this will only exist if triRep is a struct, not TriRep or Triangulation object
        data['surface']['triRep']['Triangulation']
        
    except ValueError as e:
        
        if str(e) != 'no field of name Triangulation':
            raise e
        else:
            message = 'MATLAB classes cannot be read. Please use OpenEP MATLAB to save the triRep data as a struct.'
            raise TypeError(message) from e

    # rfindex is a matlab class - not readable with Python
    data.pop('rfindex', None)

    return data
