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
    'interpolate_voltage_onto_surface',
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
    bipolar=True,
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
        bipolar (bool): If True, the bipolar traces will be returned. If False, the unipolar
            traces will be returned.
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

    electrograms = case.electric.bipolar_egm.egm if bipolar else case.electric.unipolar_egm.egm
    names = case.electric.internal_names
    local_activation_time = case.electric.annotations.local_activation_time

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
    """

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


def calculate_voltage_from_electrograms(case, buffer=50):
    """
    Calculates the peak-to-peak bipolar voltage from electrograms.

    For each mapping point, the bipolar voltage will be calculated as the
    amplitude of its corresponding electrogram during the window of interest.

    Args:
        case (Case): openep case object
        buffer (float): Amplitudes will be calculated using the window of interest
            plus/minus this buffer time.

    Returns:
        voltages (ndarray): Bipolar voltages
    """

    electrograms = case.electric.bipolar_egm.egm.copy()

    woi_times = get_woi_times(case, buffer=buffer, relative=False)
    electrograms = electrograms[:, woi_times]

    amplitudes = np.ptp(electrograms, axis=1)

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
        method (callable): method to use for interpolation. The default is
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


def interpolate_voltage_onto_surface(
        case,
        buffer=50,
        method=scipy.interpolate.RBFInterpolator,
        method_kws=None,
        max_distance=None,
):
    """Interpolate bipolar voltage onto the points of a mesh.

    For each mapping point within the window of interest, the bipolar voltage is
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

    Returns:
        interpolated_voltages (ndarray): bipolar voltages, calculated from the
        electrograms, interpolated onto the surface of the mesh.
    """

    surface_points = case.points.copy()
    points = case.electric.bipolar_egm.points.copy()
    bipolar_voltages = calculate_voltage_from_electrograms(case, buffer=buffer)

    within_woi = get_mapping_points_within_woi(case, buffer=buffer)
    interpolator = Interpolator(
        points[within_woi],
        bipolar_voltages[within_woi],
        method=method,
        method_kws=method_kws,
    )

    interpolated_voltages = interpolator(surface_points, max_distance=max_distance)

    # Any points that are not part of the mesh faces should have bipolar voltage set to NaN
    n_surface_points = surface_points.shape[0]
    not_on_surface = ~np.in1d(np.arange(n_surface_points), case.indices)
    interpolated_voltages[not_on_surface] = np.NaN

    return interpolated_voltages
