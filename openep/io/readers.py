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
    * data exported from the MATLAB implementation of OpenEP. Further, the `rfindex` data is not yet loaded.
    * data from an openCARP simulation
    * data from a Circle CVI workspace and associated dicoms


Example of loading a dataset
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Case data can be loaded as follows:

.. code:: python

    import openep
    from openep._datasets.openep_datasets import DATASET_2

    case = openep.load_openep_mat(DATASET_2)

This will load the dataset and store the information in a `Case` object.
See :class:`openep.data_structures.case.Case` for information on the attributes
and methods of `Case`.

.. autofunction:: load_openep_mat
.. autofunction:: load_opencarp

"""

import os
import scipy.io

import numpy as np

from . import _circle_cvi
from .matlab import _load_mat_v73, _load_mat_below_v73
from ..data_structures.surface import extract_surface_data, empty_fields
from ..data_structures.electric import extract_electric_data, empty_electric
from ..data_structures.ablation import extract_ablation_data, empty_ablation
from ..data_structures.case import Case

__all__ = ["load_openep_mat", "_load_mat", "load_opencarp", "load_circle_cvi"]


def _check_mat_version_73(filename):
    """Check if a MATLAB file is of version 7.3"""

    byte_stream, file_opened = scipy.io.matlab._mio._open_file(filename, appendmat=False)
    major_version, minor_version = scipy.io.matlab._miobase.get_matfile_version(byte_stream)

    return major_version == 2


def _load_mat(filename):
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
    data = _load_mat(filename)

    if name is None:
        name = os.path.basename(filename)

    points, indices, fields = extract_surface_data(data['surface'])
    electric = extract_electric_data(data['electric'])
    ablation = extract_ablation_data(data['rf']) if 'rf' in data else None

    if 'notes' in data:
        notes = np.asarray([data['notes']])[:, np.newaxis] if isinstance(data['notes'], str) else np.asarray(data['notes']).reshape(-1, 1)
    else:
        notes = np.asarray([""], dtype=str)[:, np.newaxis]

    return Case(name, points, indices, fields, electric, ablation, notes)


def load_opencarp(
    points,
    indices,
    name=None,
    fill_fields=False,
    scale_points=1,
):
    """
    Load data from an OpenCARP simulation.

    Args:
        points (str): Path to the openCARP points file.
        indices (str): Path to the openCARP element file. Currently, only triangular meshes are
            supported.
        name (str, optional): Name of the dataset. If None, the basename of the points file
            will be used as the name.
        fill_fields (bool, optional): If True, will create scalar fields filled with NaN values
            (one value per point).
        scale_points (float, optional): Scale the point positions by this number. Useful to scaling
            the units to be in mm rather than micrometre.

    Returns:
        case (Case): an OpenEP Case object that contains the points and indices.

    Note
    ----
    All other attributes of the Case object will be empty numpy arrays.

    """

    name = os.path.basename(points) if name is None else name

    points_data = np.loadtxt(points, skiprows=1)
    points_data *= scale_points
    indices_data = np.loadtxt(indices, skiprows=1, usecols=[1, 2, 3], dtype=int)  # ignore the tag for now

    size = size=points_data.size // 3 if fill_fields else 0
    fields = empty_fields(size)
    electric = empty_electric()
    ablation = empty_ablation()
    notes = np.asarray([])

    return Case(name, points_data, indices_data, fields, electric, ablation, notes)


def load_circle_cvi(filename, dicoms_directory, extract_epi=True, extract_endo=True, return_dicoms_data=False):
    """Create a pyvista.PolyData dataset from a Circle CVI workspace and stack of dicoms.

    Args:
        filename (str or pathlib.Path): Circle CVI workspace filename (.cvi42wsx).
        dicoms_directory (str or pathlib.Path): Directory containing dicoms associated with the workspace.
        extract_epi (bool, optional): Create a mesh of the epicardial surface. Default is True.
        extract_endo (bool, optional): Create a mesh of the endocardial surface. Default is True.
        return_dicoms_data (bool, optional): Whether to return a pd.DataFrame of data associated with each
            dicom. Default is False.

    Returns:
        epi_mesh (pyvista.PolyData): A mesh of the epicardium generated from the workspace file and dicoms.
        endo_mesh (pyvista.PolyData): A mesh of the endocardium generated from the workspace file and dicoms.
        dicoms (pd.DataFrame): DataFrame containing information about each dicom used to construct the mesh.
    """

    dicoms = _circle_cvi.load_dicoms(dicoms_directory=dicoms_directory)
    contour_nodes = _circle_cvi.get_contour_nodes(filename=filename)

    contour_data, dicoms_data = _circle_cvi.get_contours(
        contour_nodes,
        dicoms,
        extract_epi=extract_epi,
        extract_endo=extract_endo,
    )

    if extract_epi:

        epi_contours = [contour[key].astype(float) for contour in contour_data for key in contour if 'epi' in key ]
        for contour_index, (xy_resolution, upsample_factor) in enumerate(zip(dicoms_data.pixel_spacing_x, dicoms_data.upsample_factor)):
            epi_contours[contour_index] *= xy_resolution / upsample_factor

        epi_mesh = _circle_cvi.create_mesh(
            dicoms=dicoms_data,
            contours_xy=epi_contours,
            align_contours=True,
            n_apex_slices=2,
        )

    if extract_endo:

        endo_contours = [contour[key].astype(float) for contour in contour_data for key in contour if 'endo' in key ]
        for contour_index, (xy_resolution, upsample_factor) in enumerate(zip(dicoms_data.pixel_spacing_x, dicoms_data.upsample_factor)):
            endo_contours[contour_index] *= xy_resolution / upsample_factor

        endo_mesh = _circle_cvi.create_mesh(
            dicoms=dicoms_data,
            contours_xy=endo_contours,
            align_contours=True,
            n_apex_slices=1,
        )

    if return_dicoms_data:

        if extract_epi and extract_endo:
            return epi_mesh, endo_mesh, dicoms_data
        elif extract_epi:
            return epi_mesh, dicoms_data
        elif extract_endo:
            return endo_mesh, dicoms_data

    if extract_epi and extract_endo:
        return epi_mesh, endo_mesh
    elif extract_epi:
        return epi_mesh
    elif extract_endo:
        return endo_mesh
