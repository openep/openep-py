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
import re
import scipy.io

import numpy as np
import pyvista

from . import _circle_cvi
from .matlab import _load_mat_v73, _load_mat_below_v73
from ..data_structures.surface import extract_surface_data, Fields
from ..data_structures.electric import extract_electric_data, Electric
from ..data_structures.ablation import extract_ablation_data, Ablation
from ..data_structures.arrows import Arrows
from ..data_structures.case import Case

__all__ = ["load_openep_mat", "_load_mat", "load_opencarp", "load_circle_cvi", "load_vtk", "load_igb"]


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
    fibres=None,
    name=None,
    scale_points=1,
):
    """
    Load data from an OpenCARP simulation.

    Args:
        points (str): Path to the openCARP points file.
        indices (str): Path to the openCARP element file. Currently, only triangular meshes are
            supported.
        fibres (str, optional): Path to the openCARP fibres file.
        name (str, optional): Name of the dataset. If None, the basename of the points file
            will be used as the name.
        scale_points (float, optional): Scale the point positions by this number. Useful to scaling
            the units to be in mm rather than micrometre.

    Returns:
        case (Case): an OpenEP Case object that contains the points, indices and fibres.

    Note
    ----
    All other attributes of the Case object will be set to None.

    """

    name = os.path.basename(points) if name is None else name

    points_data = np.loadtxt(points, skiprows=1)
    points_data *= scale_points
    fibres_data = None if fibres is None else np.loadtxt(fibres)

    indices_data, cell_region_data = [], []
    linear_connection_data, linear_connection_regions = [], []

    with open(indices) as elem_file:
        data = elem_file.readlines()
        for elem in data:
            parts = elem.strip().split()
            if parts[0] == 'Tr':
                indices_data.append(list(map(int, parts[1:4])))
                cell_region_data.append(int(parts[4]))
            elif parts[0] == 'Ln':
                linear_connection_data.append(list(map(int, parts[1:3])))
                linear_connection_regions.append(int(parts[3]))

    indices_data = np.array(indices_data)
    cell_region = np.array(cell_region_data)
    linear_connection_data = np.array(linear_connection_data)
    linear_connection_regions = np.array(linear_connection_regions)

    arrows = Arrows(
        fibres=fibres_data,
        linear_connections=linear_connection_data,
        linear_connection_regions=linear_connection_regions
    )

    fields = Fields(
        cell_region=cell_region,
        longitudinal_fibres=None,
        transverse_fibres=None,
    )
    electric = Electric()
    ablation = Ablation()
    notes = np.asarray([], dtype=object)

    return Case(name, points_data, indices_data, fields, electric, ablation, notes, arrows)


def load_vtk(filename, name=None):
    """
    Load data from a VTK file.

    Args:
        filename (str): path to VTK file to be loaded (including the .vtk
            extension.)
        name (str): name to give this dataset. The default is `None`, in which case
            the filename is used at the name.

    Returns:
        case (Case): an OpenEP Case object that contains the surface, electric and
            ablation data.
    """

    name = name if name is not None else os.path.basename(filename)
    mesh = pyvista.read(filename)

    case = Case(
        name=name,
        points=mesh.points,
        indices=np.array(mesh.faces).reshape(mesh.n_cells, 4)[:, 1:],
        fields=Fields.from_pyvista(mesh),
        electric = Electric(),
        ablation = Ablation(),
        notes = np.asarray([], dtype=object),
    )

    return case


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


def load_igb(igb_filepath):
    """
    Reads an .igb file, returning the data and header information.

    Args:
        igb_filepath (str): Path to the .igb file.

    Returns:
        tuple:
            - numpy.ndarray: 2D array of the file's data.
            - dict: Contents of the header including 't' value (time steps) and other parameters.
    """
    with open(igb_filepath, 'rb') as file:
        header = file.read(1024).decode('utf-8')
        header = header.replace('\r', ' ').replace('\n', ' ').replace('\0', ' ')
        hdr_content = {}

        # Parse the header to dict format
        for part in header.split():
            key, value = part.split(':')
            if key in ['x', 'y', 'z', 't', 'bin', 'num', 'lut', 'comp']:
                hdr_content[key] = int(value)
            elif key in ['facteur','zero','epais'] or key.startswith('org_') or key.startswith('dim_') or key.startswith('inc_'):
                hdr_content[key] = float(value)
            else:
                hdr_content[key] = value

        # Process file data
        words = header.split()
        word = [int(re.split(r"(\d+)", w)[1]) for w in words[:4]]
        nnode = word[0] * word[1] * word[2]
        size = os.path.getsize(igb_filepath) // 4 // nnode

        file.seek(1024)
        data = np.fromfile(file, dtype=np.float32, count=size * nnode)
        data = data.reshape((size, nnode)).transpose()

    return data, hdr_content
