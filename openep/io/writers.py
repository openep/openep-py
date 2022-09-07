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
the points data will be stored in `dataset_2.pts`, the vertex data will be
stored in `dataset_2.elem`, and fibre orientation will be stored in `dataset_2.lon`.

Warning
-------
If a case has no fibre orientations (in `case.fields.longitudinal_fibres` and
`case.fields.transverse_fibres`), then isotropic firbres will be used.

.. autofunction:: export_openCARP

"""

import pathlib
import numpy as np
import scipy.io

from openep.data_structures.ablation import Ablation
from openep.data_structures.case import Case
from openep.data_structures.surface import Fields
from openep.data_structures.electric import Electric

__all__ = [
    "export_openCARP",
    "export_openep_mat",
    "export_vtk",
]


def export_openCARP(
    case: Case,
    prefix: str,
    scale_points: float = 1,
    export_transverse_fibres: bool = True,
):
    """Export mesh data from an OpenEP data to openCARP format.

    The following data are written to file:
        - `case.points` is used to write f'{prefix}'.pts
        - `case.indices` is used to write f'{prefix}'.elem
        - `case.fields.cell_region` is used to define the regions in f'{prefix}'.elem
        - `case.fields.longitudinal_fibres` and `case.fields.transverse_fibres` are
          used to write f'{prefix}'.lon
        - `case.fields.pacing_site` is used to write f'{pacing_site}'.vtx

    Args:
        case (Case): dataset to be exported
        prefix (str): filename prefix for writing mesh data to files
        scale_points (float, optional): Scale the point positions by this number. Useful to scaling
            the units to be in micrometre rather than mm (by setting scale_points=1e-3).
        export_transverse_fibres (bool, optional). If True, both longitudinal and transverse
            fibres directions will be written to the .lon file. If False, only longitudinal fibre
            directions will be written.
    """

    output_path = pathlib.Path(prefix).resolve()

    # Save points info
    np.savetxt(
        output_path.with_suffix(".pts"),
        case.points * scale_points,
        fmt="%.6f",
        header=str(case.points.size//3),
        comments='',
    )

    # Save elements info
    n_triangles = case.indices.shape[0]
    cell_type = np.full(n_triangles, fill_value="Tr")
    cell_region = case.fields.cell_region if case.fields.cell_region is not None else np.zeros(n_triangles, dtype=int)

    elements = np.concatenate(
        [
            cell_type[:, np.newaxis],
            case.indices,
            cell_region[:, np.newaxis],
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

    # Save fibres
    n_fibre_vectors = 2 if export_transverse_fibres else 1
    fibres = np.zeros((n_triangles, n_fibre_vectors * 3), dtype=float)

    if case.fields.longitudinal_fibres is not None:
        fibres[:, :3] = case.fields.longitudinal_fibres
    else:
        fibres[:, 0] = 1

    if export_transverse_fibres:
        if case.fields.transverse_fibres is not None:
            fibres[:, 3:] = case.fields.transverse_fibres
        else:
            fibres[:, 3] = 1

    np.savetxt(
        output_path.with_suffix('.lon'),
        fibres,
        fmt="%.6f",
        header=str(n_fibre_vectors),
        comments='',
    )

    # Saving pacing sites if they exist
    if case.fields.pacing_site is None:
        return

    # Ignore all -1 values, as these points do not belong to any pacing site
    for site_index in np.unique(case.fields.pacing_site)[1:]:
        
        pacing_site_points = np.nonzero(case.fields.pacing_site == site_index)[0]
        n_points = pacing_site_points.size

        np.savetxt(
            output_path.with_name(f'{output_path.name}_pacing_site_{site_index}.vtx'),
            pacing_site_points,
            header=f'{n_points}\nintra',
            comments='',
            fmt='%d',
        )


def export_openep_mat(
    case: Case,
    filename: str,
    separate_regions: bool = False,
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


def export_vtk(
    case: Case,
    filename: str,
):
    """Export data in vtk format.

    Saves all field data.

    Args:
        case (Case): dataset to be exported
        filename (str): name of file to be written
    """

    mesh = case.create_mesh()
    for field in case.fields:
        if case.fields[field] is None:
            continue
        if len(case.fields[field]) == mesh.n_points:
            mesh.point_data[field] = case.fields[field]
        elif len(case.fields[field]) == mesh.n_cells:
            mesh.cell_data[field] = case.fields[field]

    mesh.save(filename)


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
    surface_data['cell_region'] = fields.cell_region if fields.cell_region is not None else empty_int_array
    surface_data['fibres'] = {}
    surface_data['longitudinal'] = fields.longitudinal_fibres if fields.longitudinal_fibres is not None else empty_float_array
    surface_data['transverse'] = fields.transverse_fibres if fields.transverse_fibres is not None else empty_float_array
    surface_data['pacing_site'] = fields.pacing_site if fields.pacing_site is not None else empty_int_array

    # Remove arrays that are full of NaNs
    for field_name, field in surface_data.items():
        if isinstance(field, dict):
            continue
        if np.all(np.isnan(field)):
            surface_data[field_name] = empty_float_array

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
    electric_data['tags'] = electric._names.astype(object) if electric._names is not None else empty_object_array
    electric_data['names'] = electric._internal_names.astype(object) if electric._internal_names is not None else empty_object_array
    electric_data['include'] = electric._include if electric._include is not None else empty_int_array

    electric_data['sampleFrequency'] = float(electric.frequency)

    electric_data['electrodeNames_bip'] = electric.bipolar_egm._names.astype(object) if electric.bipolar_egm._names is not None else empty_object_array
    electric_data['egmX'] = electric.bipolar_egm._points if electric.bipolar_egm._points is not None else empty_float_array
    electric_data['egm'] = electric.bipolar_egm._egm if electric.bipolar_egm._egm is not None else empty_float_array
    electric_data['egmGain'] = electric.bipolar_egm._gain if electric.bipolar_egm._gain is not None else empty_float_array

    electric_data['electrodeNames_uni'] = electric.unipolar_egm._names.astype(object) if electric.unipolar_egm._names is not None else empty_object_array
    electric_data['egmUniX'] = electric.unipolar_egm._points if electric.unipolar_egm._points is not None else empty_float_array
    electric_data['egmUni'] = electric.unipolar_egm._egm if electric.unipolar_egm._egm is not None else empty_float_array
    electric_data['egmUniGain'] = electric.unipolar_egm._gain if electric.unipolar_egm._gain is not None else empty_float_array

    electric_data['egmRef'] = electric.reference_egm._egm if electric.reference_egm._egm is not None else empty_float_array
    electric_data['egmRefGain'] = electric.reference_egm._gain if electric.reference_egm._gain is not None else empty_float_array

    electric_data['ecg'] = electric.ecg._ecg if electric.ecg._ecg is not None else empty_float_array
    electric_data['ecgGain'] = electric.ecg._gain if electric.ecg._gain is not None else empty_float_array
    electric_data['ecgNames'] = electric.ecg.channel_names.astype(object) if electric.ecg.channel_names is not None else empty_object_array

    electric_data['egmSurfX'] = electric.surface._nearest_point if electric.surface._nearest_point is not None else empty_float_array
    electric_data['barDirection'] = electric.surface._normals if electric.surface._normals is not None else empty_float_array

    annotations = electric.annotations
    electric_data['annotations'] = {}
    electric_data['annotations']['woi'] = annotations._window_of_interest_indices if annotations._window_of_interest_indices is not None else empty_int_array
    electric_data['annotations']['referenceAnnot'] = annotations._reference_activation_time_indices if annotations._reference_activation_time_indices is not None else empty_int_array
    electric_data['annotations']['mapAnnot'] = annotations._local_activation_time_indices if annotations._local_activation_time_indices is not None else empty_int_array

    electric_data['voltages'] = {}
    electric_data['voltages']['bipolar'] = electric.bipolar_egm._voltage if electric.bipolar_egm._voltage is not None else empty_float_array
    electric_data['voltages']['unipolar'] = electric.unipolar_egm._voltage if electric.unipolar_egm._voltage is not None else empty_float_array

    # Voltages are added when loading a dataset if egms are present
    # These should be removed before saving
    for channel, voltage in electric_data['voltages'].items():
        if all(np.isnan(voltage)):
            electric_data['voltages'][channel] = empty_float_array

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