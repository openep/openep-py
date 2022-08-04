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
Saving datasets --- :mod:`openep.io.writers`
============================================

This module contains functions to export OpenEP data sets.

Example of saving a dataset
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Mesh data can be exported to openCARP format as follows:

.. code:: python

    import openep
    from openep._datasets.openep_datasets import DATASET_2

    case = openep.load_openep_mat(DATASET_2)
    openep.export_openCARP(
        case=case,
        prefix="dataset_2",
    )

This will save the mesh data from the case dataset into openCARP files, i.e.
the points data will be stored in `dataset_2.pts` and the vertex data will be
stored in `dataset_2.elem`.

Warning
-------
Currently, the fibre orientation is not calculated and written to file. Nor is the mesh
processed to make it suitable for performing openCARP simulations.

.. autofunction:: export_openCARP

"""

import pathlib
from idna import alabel
import numpy as np
import scipy.io
from pyvista import PolyData

from openep.data_structures.ablation import Ablation
from openep.data_structures.case import Case
from openep.data_structures.surface import Fields
from openep.data_structures.electric import Electric

__all__ = ["export_openCARP", "export_openep_mat"]


def export_openCARP(
    case: Case,
    prefix: str,
):
    """Export mesh data from an OpenEP data to openCARP format.

    Args:
        case (Case): dataset to be exported
        prefix (str): filename prefix for writing mesh data to files
    """

    output_path = pathlib.Path(prefix).resolve()

    # Save points info
    np.savetxt(
        output_path.with_suffix(".pts"),
        case.points,
        fmt="%.6f",
        header=str(case.points.size),
        comments='',
    )

    # Save elements info
    n_triangles = case.indices.shape[0]
    cell_type = np.full(n_triangles, fill_value="Tr")
    region = np.full(n_triangles, 100, dtype=int)

    elements = np.concatenate(
        [
            cell_type[:, np.newaxis],
            case.indices,
            region[:, np.newaxis],
        ],
        axis=1,
        dtype=object,
    )

    np.savetxt(
        output_path.with_suffix(".elem"),
        elements,
        fmt="%s %d %d %d %d",
        header=str(n_triangles),
        comments='',
    )


def export_openep_mat(
    case: Case,
    filename: str,
):
    """Export data in OpenEP format.

    Args:
        case (Case): dataset to be exported
        filename (str): name of file to be written
    """

    userdata = {}
    userdata['notes'] = case.notes.astype(object) if case.notes is not None else np.array([], dtype=object)

    userdata['surface'] = _extract_surface_data(
        points=case.points,
        indices=case.indices,
        fields=case.fields,
    )

    userdata['electric'] = _extract_electric_data(electric=case.electric)
    userdata['rf'] = _export_ablation_data(ablation=case.ablation)

    scipy.io.savemat(
        file_name=filename,
        mdict={'userdata': userdata},
        format="5",
        do_compression=True,
        oned_as='column',
    )

def _extract_surface_data(
    points : np.ndarray,
    indices: np.ndarray,
    fields: Fields,
):
    """Create a dictionary of surface data.

    Args:
        points (np.ndarray): 3D coordinates of points on the mesh
        indices (np.ndarray): Indices of points that make up each face of the mesh
        fields (Fields): bipolar voltage, unipolar voltage, local activation time,
            impedance, and force at each point on the mesh.

    Returns:            
        surface_data (dict): Dictionary containing numpy arrays that describe the
            surface of a mesh as well as scalar values (fields)
    """

    # Make sure none of the attributes are equal to None - scipy can't handle this
    empty_float_array = np.array([], dtype=float)
    empty_int_array = np.array([], dtype=int)

    surface_data = {}

    surface_data['triRep'] = {}
    surface_data['triRep']['X'] = points if points is not None else empty_float_array

    # TriRep requires doubles not ints, and MATLAB uses 1-based indexing
    triangulation = indices.astype(float) + 1 if indices is not None else empty_int_array
    surface_data['triRep']['Triangulation'] = triangulation

    if fields.local_activation_time is None and fields.bipolar_voltage is None:
        surface_data['act_bip'] = empty_float_array
    elif fields.local_activation_time is None:
        fields.local_activation_time = np.full_like(fields.bipolar_voltage, fill_value=np.NaN)
    elif fields.bipolar_voltage is None:
        fields.bipolar_voltage = np.full_like(fields.local_activation_time, fill_value=np.NaN)

    if 'act_bip' not in surface_data:
        surface_data['act_bip'] = np.concatenate(
            [
                fields.local_activation_time[:, np.newaxis],
                fields.bipolar_voltage[:, np.newaxis],
            ],
            axis=1,
        )

    if fields.unipolar_voltage is None and fields.impedance is None and fields.force is None:
        surface_data['uni_imp_frc'] = empty_float_array
    else:
        if fields.unipolar_voltage is None:
            fields.unipolar_voltage = np.full(points.size // 3, fill_value=np.NaN)
        if fields.impedance is None:
            fields.impedance = np.full(points.size // 3, fill_value=np.NaN)
        if fields.force is None:
            fields.force = np.full(points.size // 3, fill_value=np.NaN)

    if 'uni_imp_frc' not in surface_data:
        surface_data['uni_imp_frc'] = np.concatenate(
            [
                fields.unipolar_voltage[:, np.newaxis],
                fields.impedance[:, np.newaxis],
                fields.force[:, np.newaxis],
            ],
            axis=1,
        )

    surface_data['thickness'] = fields.thickness if fields.thickness is not None else empty_float_array

    return surface_data


def _extract_electric_data(electric: Electric):
    """Create a dictionary of electric data.

    Args:
        electric (Electric): object containing electric data associated with electrograms
            taken at various mapping points.

    Returns:
        electric_data (dict): Dictionary containing numpy arrays that describe the
            electric data associated with electrograms taken at various mapping points.
    """

    # Make sure none of the attributes are equal to None - scipy can't handle this
    empty_object_array = np.array([], dtype=object)  # for empty string arrays
    empty_float_array = np.array([], dtype=float)
    empty_int_array = np.array([], dtype=int)

    electric_data = {}
    electric_data['tags'] = electric.names.astype(object) if electric.names is not None else empty_object_array
    electric_data['names'] = electric.internal_names.astype(object) if electric.internal_names is not None else empty_object_array
    electric_data['include'] = electric.include if electric.include is not None else empty_int_array

    electric_data['sampleFrequencu'] = float(electric.frequency)

    electric_data['electrodeNames_bip'] = electric.bipolar_egm.names.astype(object) if electric.bipolar_egm.names is not None else empty_object_array
    electric_data['egmX'] = electric.bipolar_egm.points if electric.bipolar_egm.points is not None else empty_float_array
    electric_data['egm'] = electric.bipolar_egm.egm if electric.bipolar_egm.egm is not None else empty_float_array
    electric_data['egmGain'] = electric.bipolar_egm.gain if electric.bipolar_egm.gain is not None else empty_float_array

    electric_data['electrodeNames_uni'] = electric.unipolar_egm.names.astype(object) if electric.unipolar_egm.names is not None else empty_object_array
    electric_data['egmUniX'] = electric.unipolar_egm.points if electric.unipolar_egm.points is not None else empty_float_array
    electric_data['egmUni'] = electric.unipolar_egm.egm if electric.unipolar_egm.egm is not None else empty_float_array
    electric_data['egmUniGain'] = electric.unipolar_egm.gain if electric.unipolar_egm.gain is not None else empty_float_array

    electric_data['egmRef'] = electric.reference_egm.egm if electric.reference_egm.egm is not None else empty_float_array
    electric_data['egmRefGain'] = electric.reference_egm.gain if electric.reference_egm.gain is not None else empty_float_array

    electric_data['ecg'] = electric.ecg.ecg if electric.ecg.ecg is not None else empty_float_array
    electric_data['ecgGain'] = electric.ecg.gain if electric.ecg.gain is not None else empty_float_array
    electric_data['ecgNames'] = electric.ecg.channel_names.astype(object) if electric.ecg.channel_names is not None else empty_object_array

    electric_data['egmSurfX'] = electric.surface.nearest_point if electric.surface.nearest_point is not None else empty_float_array
    electric_data['barDirection'] = electric.surface.normals if electric.surface.normals is not None else empty_float_array

    electric_data['annotations'] = {}
    electric_data['annotations']['woi'] = electric.annotations.window_of_interest
    electric_data['annotations']['referenceAnnot'] = electric.annotations.reference_activation_time
    electric_data['annotations']['mapAnnot'] = electric.annotations.local_activation_time

    electric_data['voltages'] = {}
    electric_data['voltages']['bipolar'] = electric.bipolar_egm.voltage if electric.bipolar_egm.voltage is not None else empty_float_array
    electric_data['voltages']['unipolar'] = electric.unipolar_egm.voltage if electric.unipolar_egm.voltage is not None else empty_float_array

    electric_data['impedances'] = {}
    electric_data['impedances']['time'] = electric.impedance.times if electric.impedance.times is not None else empty_float_array
    electric_data['impedances']['value'] = electric.impedance.values if electric.impedance.values is not None else empty_float_array

    return electric_data


def _export_ablation_data(ablation: Ablation):
    """Create a dictionary of ablation data.

    Args:
        ablation (Ablation): times, power, impedance and temperature for each ablation site,
            as well as the force applied.

    Returns:
        ablation_data (dict): Dictionary containing numpy arrays that describe the
            ablation sites.
    """

    # Make sure none of the attributes are equal to None - scipy can't handle this
    empty_float_array = np.array([], dtype=float)

    ablation_data = {}
    ablation_data['originaldata'] = {}

    ablation_data['originaldata']['ablparams'] = {}
    ablation_data['originaldata']['ablparams']['time'] = ablation.times if ablation.times is not None else empty_float_array
    ablation_data['originaldata']['ablparams']['power'] = ablation.power if ablation.power is not None else empty_float_array
    ablation_data['originaldata']['ablparams']['impedance'] = ablation.impedance if ablation.impedance is not None else empty_float_array
    ablation_data['originaldata']['ablparams']['distaltemp'] = ablation.temperature if ablation.temperature is not None else empty_float_array

    ablation_data['originaldata']['force'] = {}
    ablation_data['originaldata']['force']['time'] = ablation.force.times if ablation.force.times is not None else empty_float_array
    ablation_data['originaldata']['force']['force'] = ablation.force.force if ablation.force.force is not None else empty_float_array
    ablation_data['originaldata']['force']['axialangle'] = ablation.force.axial_angle if ablation.force.axial_angle is not None else empty_float_array
    ablation_data['originaldata']['force']['lateralangle'] = ablation.force.lateral_angle if ablation.force.lateral_angle is not None else empty_float_array
    ablation_data['originaldata']['force']['position'] = ablation.force.points if ablation.force.points is not None else empty_float_array

    return ablation_data