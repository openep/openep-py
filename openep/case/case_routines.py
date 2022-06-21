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
Analyse a Case dataset - :mod:`openep.case.case_routines`
=========================================================

This module provides methods for :ref:`analysing <analysing>` or
:ref:`extracting data <extracting>` from a case object, as
well as :ref:`interpolating <interpolating>` electrical data from
mapping points onto the 3D surface.


.. _extracting:

Extracting data from a Case
---------------------------

.. autofunction:: get_electrograms_at_points

.. autofunction:: get_woi_times

.. _analysing:

Analyses
--------

.. autofunction:: calculate_voltage_from_electrograms

.. _interpolating:

Interpolate electrical data from mapping points onto the 3D surface
-------------------------------------------------------------------


.. autofunction:: interpolate_voltage_onto_surface

Tip
---

    By default, `scipy`'s `RBFInterpolator` is used to perform the interpolation.
    If you prefer to use a different interpolation method, you pass an interpolator
    class and associated keyword arguments to the `method` and `method_kws`
    arugments respectively.

.. autoclass:: Interpolator
    :members: __call__

"""

from attr import attrs

import numpy as np
import scipy.interpolate

__all__ = [
    'get_mapping_points_within_woi',
    'get_electrograms_at_points',
    'get_woi_times',
    'calculate_voltage_from_electrograms',
    'calculate_distance',
    'calculate_points_within_distance',
    'Interpolator',
    'interpolate_activation_time_onto_surface',
    'interpolate_voltage_onto_surface',
    'bipolar_from_unipolar_surface_points',
]


def _get_reference_annotation(case, indices=None):
    """
    Get the reference activation time for each mapping point.

    Args:
        case (Case): openep case object
        indices (ndarray), optional: indices of reference annotations to select. The default
            is None, in which case all reference annotation will be returned.
    Returns:
        annotations (ndarray): reference annotations
    """

    annotations = case.electric.annotations.reference_activation_time
    annotations = annotations[indices] if indices is not None else annotations

    return annotations


def _get_window_of_interest(case, indices=None):
    """
    Gets the window of interest for each mapping point.

    Args:
        case (Case): openep case object
        indices (ndarray), optional: indices of mapping points for which the woi will
            be extracted. The default is None, in which case the woi will be returned for
            every mapping point.

    Returns:
        woi (ndarray): Nx2 array with the windows of interest for selected mapping points.
            For each point, the two columns specify that start and end time respectively
            of the window of interest.
    """

    woi = case.electric.annotations.window_of_interest
    woi = woi[indices] if indices is not None else woi.copy()

    return woi


def get_mapping_points_within_woi(case, indices=None, buffer=50):
    """
    Extracts the indices of the mapping points that have
    annotated local activation times within the window of interest (woi).

    Args:
        case (Case): openep case object
        indices (ndarray), optional: indices of mapping points to consider. The default
            is None, in which case all points within the window of interest will be returned.
        buffer (float): points within the woi plus/minus this buffer time will
            be considered to be within the woi.

    Returns:
        within_woi (ndarray): boolean array of valid points; True if
            the annotated local activation time is within the woi, False otherwise.
    """

    # if we have a single index we need to ensure it is an array
    indices = np.asarray([indices], dtype=int) if isinstance(indices, int) else indices

    reference_activation_times = _get_reference_annotation(case, indices=indices)
    woi = _get_window_of_interest(case, indices=indices)
    woi += reference_activation_times[:, np.newaxis]

    local_activation_times = case.electric.annotations.local_activation_time
    local_activation_times = local_activation_times[indices] if indices is not None else local_activation_times

    within_woi = np.logical_and(
        local_activation_times >= woi[:, 0] - buffer,
        local_activation_times <= woi[:, 1] + buffer,
    )

    return within_woi


def get_electrograms_at_points(
    case,
    within_woi=True,
    buffer=50,
    indices=None,
    egm_type="bipolar",
    return_names=True,
    return_lat=True,
):
    """
    Extract the electrogram timeseries for a selection of points.

    Args:
        case (Case): openep case object
        within_woi (bool): If True, only electrograms within the window of interest will be extracted.
            If False, all electrograms will be extracted.
        buffer (float): If `within_woi` is True, the woi will be extended by this time. If `wihtin_woi`
            is False, `buffer` is ignored.
        indices (ndarray), optional: indices of mapping points for which electrograms will
            be extracted. If provided along with `woi=True`, only the electrograms that
            are both within the window of interest and selected by `indices` will be extracted.
            If provided along with `woi=False`, all electrograms of mapping points selected by
            `indices` will be returned.
        egm_type (bool): The electrogram traces to return. Valid options:
            - bipolar
            - unipolar
            - reference
        return_names (bool): If True, the internal names of the points used by the
            clinical electroanatomic mapping system will also be returned.
        return_lat (bool): If True, the local activation time of each mapping point will also be
            returned.

    Returns:
        traces (ndarray): A timeseries of voltages for each selected mapping point.

        names (ndarray, optional): If `return_names` is True, the internal names of the points used by the
        clinical electroanatomic mapping system will be returned.

        local_activation_time (ndarray, optional): If `return_lat` is True, the local activation time of
        each mapping point will be returned.

    """

    if egm_type.lower().strip() == "bipolar":
        electrograms = case.electric.bipolar_egm.egm
    elif egm_type.lower().strip() == "unipolar":
        electrograms = case.electric.unipolar_egm.egm
    elif egm_type.lower().strip() == "reference":
        electrograms = case.electric.reference_egm.egm
    else:
        raise ValueError(f"egm_type {egm_type} is not recognised.")

    if case.electric.internal_names is not None:
        names = case.electric.internal_names
    else:
        names = np.asarray([f"P{index}" for index in range(len(electrograms))], dtype=str)

    if case.electric.annotations.local_activation_time is not None:
        local_activation_time = case.electric.annotations.local_activation_time
    else:
        local_activation_time = np.full_like(names, fill_value=np.NaN)

    # Filter by selected indices
    if indices is not None:

        # if we have a single index we need to ensure it is an array
        indices = np.asarray([indices], dtype=int) if isinstance(indices, int) else indices

        electrograms = electrograms[indices]
        names = names[indices]
        local_activation_time = local_activation_time[indices]

    # Filter by window of interest and buffer
    if within_woi:
        points_within_woi = get_mapping_points_within_woi(case, indices=indices, buffer=buffer)
        electrograms = electrograms[points_within_woi]
        names = names[points_within_woi]
        local_activation_time = local_activation_time[points_within_woi]

    if return_names and return_lat:
        return electrograms, names, local_activation_time
    elif return_names:
        return electrograms, names
    elif return_lat:
        return electrograms, local_activation_time

    return electrograms


def get_woi_times(case, buffer=50, relative=False):
    """
    Extracts the times within the window of interest.

    Args:
        case (Case): openep case object
        buffer (float): times within the window of interest plus/minus this buffer
            time will be considered to be within the woi.
        relative (bool): if True, the time of the reference annotation will be
            subtracted from the returned times.

    Returns:
        times (ndarray): times within the window of interest
        
    Warning
    -------
    This returns the sample indices that cover the window of interest for the
    **first** electrogram only. get_sample_indices_within_woi can be used to
    obtain a two-dimensional boolean array in which value of True indicate
    that the specific sample is within the specific electrogram's window of
    interest.
    """
    
    # TODO: remove this function as it does not take into account the fact
    #       that different electrograms may have different windows of interest
    #       and different reference annotations

    # TODO: This is sample rather than time
    #       Should take into account sample frequency to get actual time
    times = np.arange(case.electric.bipolar_egm.egm.shape[1])

    # TODO: woi and reference times might be different for each mapping point
    woi = case.electric.annotations.window_of_interest[0]
    ref_annotation = case.electric.annotations.reference_activation_time[0]

    start_time, stop_time = woi + ref_annotation + [-buffer, buffer]

    keep_times = np.logical_and(
        times >= start_time,
        times <= stop_time,
    )

    if relative:
        times -= ref_annotation

    return times[keep_times]


def get_sample_indices_within_woi(case, buffer=50):
    """
    Determine which samples are within the window of interest for each electrogram.

    can be used to
    obtain a two-dimensional boolean array in which value of True indicate
    that the specific sample is within the specific electrogram's window of
    interest.

    Args:
        case (Case): openep case object
        buffer (float): times within the window of interest plus/minus this buffer
            time will be considered to be within the woi.

    Returns:
        within_woi (ndarray): 2D boolean array of shape (N_electrograms, N_samples).
            Values of True indicate that the sample if within the electrogram's
            window of interest.
    """
    
    egm = case.electric.bipolar_egm.egm
    sample_indices = np.full_like(egm, fill_value=np.arange(egm.shape[1]), dtype=int)
    
    woi = case.electric.annotations.window_of_interest
    ref_annotations = case.electric.annotations.reference_activation_time[:, np.newaxis]

    start_time, stop_time = (woi + ref_annotations + [-buffer, buffer]).T

    within_woi = np.logical_and(
        sample_indices >= start_time[:, np.newaxis],
        sample_indices <= stop_time[:, np.newaxis],
    )

    return within_woi  # This is now a 2D array that can be used to index into electrograms and calculate voltages.


def calculate_voltage_from_electrograms(case, buffer=50, bipolar=True):
    """
    Calculates the peak-to-peak voltage from electrograms.

    For each mapping point, the voltage will be calculated as the
    amplitude of its corresponding electrogram during the window of interest.

    Args:
        case (Case): openep case object
        buffer (float): Amplitudes will be calculated using the window of interest
            plus/minus this buffer time.
        bipolar (bool, optional): If True, the bipolar voltages will calculated. If False,
            the unipolar voltages will be calculated.

    Returns:
        voltages (ndarray): Bipolar voltages
    """

    if bipolar:
        electrograms = case.electric.bipolar_egm.egm.copy()
    else:
        electrograms = case.electric.unipolar_egm.egm.copy()

    sample_within_woi = get_sample_indices_within_woi(case, buffer=buffer)
    electrograms[~sample_within_woi] = np.NaN

    amplitudes = np.nanmax(electrograms, axis=1) - np.nanmin(electrograms, axis=1)

    return amplitudes


def calculate_distance(origin, destination):
    """
    Returns the distance from a set of origin points to a set of destination
    points.

    Args:
        origin (ndarray): Nx3 matrix of coordinates
        destination (ndarray): Mx3 matrix of coordinates

    Returns:
        distances (ndarray): MxN matrix of distances

    """

    origin = origin[np.newaxis, :] if origin.ndim == 1 else origin
    destination = destination[np.newaxis, :] if destination.ndim == 1 else destination

    distances = scipy.spatial.distance.cdist(
        origin,
        destination,
    )

    return distances


def calculate_points_within_distance(origin, destination, max_distance, return_distances=True):
    """
    Calculates whether the distances from a set of origin points to a set of
    destination points are each within a specified distance cutoff.

    Args:
        points (ndarray, (N, 3)): array of 3D points
        test_points (ndarray, (M, 3)): array of 3D test points
        max_distance (float): distance threshold between the origin points and
            destination points

    Returns:
        within_max_dist (ndarray, M x N): Boolean array
        that is equal to True if the points are within the maximum distance
        of one another and equal to False otherwise.

        distances (ndarray, M x N): distance between each point
        and each test_point.
    """

    distances = calculate_distance(origin, destination)
    within_max_distance = distances <= max_distance

    if return_distances:
        return within_max_distance, distances

    return within_max_distance


@attrs(auto_attribs=True, auto_detect=True)
class Interpolator:
    """Interpolate scalar values onto the points of a mesh.

    Args:
        points (np.ndarray): (N,3) array of coordinates for which we know values
            of the scalar field
        field (np.ndarray): array of size N of scalar values
        method (callable): method to use for interpolation. Must be a callable class
            that performs interpolation when called. The default is
            scipy.interpolate.RBFInterpolator.
        method_kws (dict): dictionary of keyword arguments to pass to `method`
            when creating the interpolator
    """

    points: np.ndarray
    field: np.ndarray
    method: callable = scipy.interpolate.RBFInterpolator
    method_kws: dict = None

    def __attrs_post_init__(self):
        """Create the interpolator"""

        default_rbf_kws = {
            "smoothing": 4,
            "kernel": 'multiquadric',
            "epsilon": 1,
            "degree": 1,
        }

        if self.method_kws is not None:
            if self.method is scipy.interpolate.RBFInterpolator:
                self.method_kws = {**default_rbf_kws, **self.method_kws}
        elif self.method is scipy.interpolate.RBFInterpolator:
            self.method_kws = default_rbf_kws
        else:
            self.method_kws = {}

        self.interpolate = self.method(
            self.points,
            self.field,
            **self.method_kws,
        )

    def __call__(self, surface_points, max_distance=None):
        """Interpolate the scalar field onto a new set of coordinates

        Args:
            surface_points (ndarray): Mx3 array of points onto which the field will be
                interpolated
            max_distance (float, optional): Surface points further than this distance from any of the
                original points will not be used in the interpolation. Instead, their scalar field will
                be set to NaN. Defaults to None, in which case all surface points will used.

        Returns:
            interpolated_field (ndarray): Scalar field interpolated onto the new points.
        """

        interpolated_field = self.interpolate(surface_points)

        if max_distance is not None:
            within_distance = calculate_points_within_distance(
                surface_points,
                self.points,
                max_distance=max_distance,
                return_distances=False
            )
            within_distance = np.any(within_distance, axis=1)
            interpolated_field[~within_distance] = np.NaN

        return interpolated_field

    def __repr__(self):
        return f"Interpolator: method={self.method}, kws={self.method_kws}"


def interpolate_activation_time_onto_surface(
        case,
        buffer=50,
        method=scipy.interpolate.RBFInterpolator,
        method_kws=None,
        max_distance=None,
):
    """Interpolate local activation times onto the points of a mesh.

    Only mapping points within the window of interest, plus the buffer time,
    are used for the interpolation.

    Args:
        case (openep.case.Case): case from which the activation time will be interpolated
        buffer (float, optional): extend the window of interest by this time.
            The default is 50 ms.
        method (callable): method to use for interpolation. The default is
            scipy.interpolate.RBFInterpolator.
        method_kws (dict): dictionary of keyword arguments to pass to `method`
            when creating the interpolator.
        max_distance (float, optional): If provided, any points on the surface of the mesh
            further than this distance to all mapping coordinates will have their
            interpolated activation times set NaN. The default it None, in which case
            the distance from surface points to mapping points is not considered.

    Returns:
        interpolated_lat (ndarray): local activation times interpolated onto the surface of the mesh,
            one value per point on the mesh.
    """

    surface_points = case.points.copy()

    points = case.electric.bipolar_egm.points.copy()
    local_activation_times = case.electric.annotations.local_activation_time - case.electric.annotations.reference_activation_time

    within_woi = get_mapping_points_within_woi(case, buffer=buffer)
    points = points[within_woi]
    local_activation_times = local_activation_times[within_woi]

    interpolator = Interpolator(
        points,
        local_activation_times,
        method=method,
        method_kws=method_kws,
    )

    interpolated_lat = interpolator(surface_points, max_distance=max_distance)

    # Any points that are not part of the mesh faces should have bipolar voltage set to NaN
    n_surface_points = surface_points.shape[0]
    not_on_surface = ~np.in1d(np.arange(n_surface_points), case.indices)
    interpolated_lat[not_on_surface] = np.NaN

    return interpolated_lat


def interpolate_voltage_onto_surface(
        case,
        buffer=50,
        method=scipy.interpolate.RBFInterpolator,
        method_kws=None,
        max_distance=None,
        bipolar=True,
):
    """Interpolate voltage onto the points of a mesh.

    For each mapping point within the window of interest, the voltage is
    calculated as the amplitude of an electrogram during the window of interest,
    plus/minus an optional buffer time.

    Args:
        case (openep.case.Case): case from which the voltage will be calculated
        buffer (float, optional): extend the window of interest by this time.
            The default is 50 ms.
        method (callable): method to use for interpolation. The default is
            scipy.interpolate.RBFInterpolator.
        method_kws (dict): dictionary of keyword arguments to pass to `method`
            when creating the interpolator.
        max_distance (float, optional): If provided, any points on the surface of the mesh
            further than this distance to all mapping coordinates will have their
            interpolated voltages set NaN. The default it None, in which case
            the distance from surface points to mapping points is not considered.
        bipolar (bool, optional): If True, the bipolar voltage will be interpolated onto the
            surface. If False, the unipolar voltage will be used instead.

    Returns:
        interpolated_voltages (ndarray): bipolar voltages, calculated from the
        electrograms, interpolated onto the surface of the mesh.
    """

    surface_points = case.points.copy()

    if bipolar:
        points = case.electric.bipolar_egm.points.copy()
        voltages = calculate_voltage_from_electrograms(case, buffer=buffer)
    else:
        points = case.electric.unipolar_egm.points.copy()
        voltages = calculate_voltage_from_electrograms(case, buffer=buffer, bipolar=False)

    within_woi = get_mapping_points_within_woi(case, buffer=buffer)
    if bipolar:
        points = points[within_woi]
        voltages = voltages[within_woi]
    else:
        points = np.concatenate(points[within_woi], axis=1).T
        voltages = voltages[within_woi].flatten()

    interpolator = Interpolator(
        points,
        voltages,
        method=method,
        method_kws=method_kws,
    )

    interpolated_voltages = interpolator(surface_points, max_distance=max_distance)

    # Any points that are not part of the mesh faces should have bipolar voltage set to NaN
    n_surface_points = surface_points.shape[0]
    not_on_surface = ~np.in1d(np.arange(n_surface_points), case.indices)
    interpolated_voltages[not_on_surface] = np.NaN

    return interpolated_voltages


def bipolar_from_unipolar_surface_points(unipolar, indices):
    """Calculate bipolar electrograms from unipolar electrograms for each point on a mesh.

    Warning
    -------
    The number of unipolar traces **must** be equal to the number of points on the surface.
    This function should therefore only be used for openCARP datasets that have had phi recovery performed
    for each point on the surface.

    Args:
        unipolar (np.ndarray): Unipolar electrograms
        indices (np.ndarray): Indices of each point in each triangle (i.e Case.indices).

    Returns
        bipolar (np.ndarray):
            Bipolar electrograms
        pair_indices (np.ndarray):
            Indices of unipolar electrograms that contribute to each bipolar electrogram.

    """

    bipolar = np.full_like(unipolar, fill_value=np.NaN)
    pair_indices = np.full((len(unipolar), 2), fill_value=0, dtype=int)

    for index, index_unipolar in enumerate(unipolar):

        connected_vertices = _find_connected_vertices(indices, index)
        index_bipolar, pair_index = _bipolar_from_unipolar(
            unipolar=index_unipolar,
            neighbours=unipolar[connected_vertices],
        )
        bipolar[index] = index_bipolar
        pair_indices[index] = [index, connected_vertices[pair_index]]

    return bipolar, pair_indices


def _find_connected_vertices(indices, index):
    """
    Find all points connected to a given point by a single edge.

    Adapted from: https://github.com/pyvista/pyvista-support/issues/96#issuecomment-571864471

    Args:
        indices (np.ndarray): triangular faces of a mesh
        index (int): index of point for which we want to find the neighbouring points
    """

    connected_faces = [i for i, face in enumerate(indices) if index in face]
    connected_vertices = np.unique(indices[connected_faces])

    return connected_vertices[connected_vertices != index]


def _bipolar_from_unipolar(unipolar, neighbours):
    """
    Calculate the bipolar electrogram of a given point from a series of unipolar electrograms.

    Args:
        unipolar (np.ndarray): unipolar electrogram at a given point
        neighbours (np.narray): unipolar electrograms at all points
            neighbouring the given point
    """

    difference = neighbours - unipolar
    voltage = np.ptp(difference, axis=1)
    pair_index = np.argmax(voltage)
    bipolar_electrogram = difference[pair_index]

    return bipolar_electrogram, pair_index
