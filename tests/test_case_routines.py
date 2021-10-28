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

import pytest
from numpy.testing import assert_allclose, assert_array_equal

import numpy as np
import scipy.interpolate

import openep
from openep.case.case_routines import (
    _get_reference_annotation,
    _get_window_of_interest,
    get_mapping_points_within_woi,
    get_electrograms_at_points,
    get_woi_times,
    calculate_voltage_from_electrograms,
    calculate_distance,
    calculate_points_within_distance,
    Interpolator,
    interpolate_voltage_onto_surface,
)
from openep._datasets.openep_datasets import DATASET_2_V73


@pytest.fixture(scope='module')
def real_case():
    return openep.load_case(DATASET_2_V73)

@pytest.fixture()
def mock_case(mocker):
    
    case = mocker.patch('openep.data_structures.case.Case')

    case.electric.annotations.reference_activation_time = np.full(10, fill_value=5, dtype=int)
    case.electric.annotations.local_activation_time = np.arange(5, 25, 2)
    case.electric.annotations.window_of_interest = np.full((10, 2), fill_value=5)
    case.electric.annotations.window_of_interest[:, 1] += 10

    case.electric.bipolar_egm.egm = np.repeat(np.arange(20)[np.newaxis, :], 10, axis=0)
    case.electric.internal_names = np.arange(10).astype(str)

    case.points = np.repeat(np.arange(20)[:, np.newaxis], repeats=3, axis=1)
    case.electric.bipolar_egm.points = np.repeat(np.arange(10)[:, np.newaxis], repeats=3, axis=1)

    return case


def test_get_reference_annotation(mock_case):

    reference_annotations = _get_reference_annotation(mock_case, indices=None)

    assert_allclose(
        mock_case.electric.annotations.reference_activation_time,
        reference_annotations,
    )


@pytest.mark.parametrize('indices', [1, [0, 1, 5]])
def test_get_reference_annotation_indices(mock_case, indices):

    reference_annotations = _get_reference_annotation(mock_case, indices=indices)

    assert_allclose(
        mock_case.electric.annotations.reference_activation_time[indices],
        reference_annotations,
    )


def test_get_window_of_interest(mock_case):

    window_of_interest = _get_window_of_interest(mock_case, indices=None)

    assert_allclose(
        mock_case.electric.annotations.window_of_interest,
        window_of_interest,
    )


@pytest.mark.parametrize('indices', [1, [0, 1, 5]])
def test_get_window_of_interest_indices(mock_case, indices):

    window_of_interest = _get_window_of_interest(mock_case, indices=indices)

    assert_allclose(
        mock_case.electric.annotations.window_of_interest[indices],
        window_of_interest,
    )


def test_get_mapping_points_within_woi(mock_case):

    within_woi = get_mapping_points_within_woi(mock_case)

    assert within_woi.all()


def test_get_mapping_points_within_woi_indices(mock_case):

    within_woi = get_mapping_points_within_woi(mock_case, indices=[0, 1, 5])

    assert within_woi.size == 3
    assert within_woi.all()


def test_get_mapping_points_within_woi_no_buffer(mock_case):

    within_woi = get_mapping_points_within_woi(mock_case, indices=[3, 4, 5, 6, 7], buffer=0)
    outside_woi = get_mapping_points_within_woi(mock_case, indices=[0, 1, 2, 8, 9], buffer=0)

    assert within_woi.size == 5
    assert within_woi.all()

    assert outside_woi.size == 5
    assert not outside_woi.any()


def test_get_electrograms_at_points(mock_case):

    electrograms, names, lat = get_electrograms_at_points(mock_case, within_woi=False)

    assert_allclose(mock_case.electric.bipolar_egm.egm, electrograms)
    assert_array_equal(mock_case.electric.internal_names, names)
    assert_allclose(mock_case.electric.annotations.local_activation_time, lat)

def test_get_electrograms_at_points_within_woi(mock_case):

    electrograms, names, lat = get_electrograms_at_points(mock_case, within_woi=True, buffer=5)

    assert_allclose(mock_case.electric.bipolar_egm.egm, electrograms)
    assert_array_equal(mock_case.electric.internal_names, names)
    assert_allclose(mock_case.electric.annotations.local_activation_time, lat)


def test_get_electrograms_at_points_within_woi_no_buffer(mock_case):

    electrograms, names, lat = get_electrograms_at_points(mock_case, within_woi=True, buffer=0)
    within_woi = [3, 4, 5, 6, 7]

    assert_allclose(mock_case.electric.bipolar_egm.egm[within_woi], electrograms)
    assert_array_equal(mock_case.electric.internal_names[within_woi], names)
    assert_allclose(mock_case.electric.annotations.local_activation_time[within_woi], lat)


def test_get_electrograms_at_points_indices(mock_case):

    indices = [0, 1, 5]
    electrograms, names, lat = get_electrograms_at_points(mock_case, within_woi=False, indices=indices)

    assert_allclose(mock_case.electric.bipolar_egm.egm[indices], electrograms)
    assert_array_equal(mock_case.electric.internal_names[indices], names)
    assert_allclose(mock_case.electric.annotations.local_activation_time[indices], lat)


def test_get_electrograms_at_points_no_lat(mock_case):

    electrograms, names = get_electrograms_at_points(mock_case, within_woi=False, return_lat=False)

    assert_allclose(mock_case.electric.bipolar_egm.egm, electrograms)
    assert_array_equal(mock_case.electric.internal_names, names)


def test_get_electrograms_at_points_no_names(mock_case):

    electrograms, lat = get_electrograms_at_points(mock_case, within_woi=False, return_names=False)

    assert_allclose(mock_case.electric.bipolar_egm.egm, electrograms)
    assert_allclose(mock_case.electric.annotations.local_activation_time, lat)


def test_get_electrograms_at_points_no_lat_or_names(mock_case):

    electrograms = get_electrograms_at_points(mock_case, within_woi=False, return_lat=False, return_names=False)
    assert_allclose(mock_case.electric.bipolar_egm.egm, electrograms)


def test_get_woi_times(mock_case):
    
    times = get_woi_times(mock_case)
    assert_allclose(np.arange(20), times)


def test_get_woi_times_no_buffer(mock_case):
    
    times = get_woi_times(mock_case, buffer=0)
    assert_allclose(np.arange(10, 20), times)

def test_get_woi_times_no_buffer_relative(mock_case):
    
    times = get_woi_times(mock_case, buffer=0, relative=True)
    assert_allclose(np.arange(5, 15), times)


def test_calculate_voltage_from_electrograms(mock_case):

    amplitudes = calculate_voltage_from_electrograms(mock_case)
    assert_allclose(19, amplitudes)


def test_calculate_voltage_from_electrograms_no_buffer(mock_case):

    amplitudes = calculate_voltage_from_electrograms(mock_case, buffer=0)
    assert_allclose(9, amplitudes)


def test_calculate_distance(mock_case):
    
    origin = mock_case.electric.bipolar_egm.points
    destination = mock_case.points

    distances = calculate_distance(origin, destination)
    diagonal_distacnces = 0  # row-wise, origin and destination are the same

    assert_allclose(diagonal_distacnces, distances.diagonal())
    assert_allclose(distances[0, :10], distances[:, 0])  # there are only 10 emg points but 20 surface points
    assert (10, 20) == distances.shape


def test_calculate_distance_single_point(mock_case):
    
    origin = mock_case.electric.bipolar_egm.points[0]
    destination = mock_case.points[-1]

    distances = calculate_distance(origin, destination)
    actual_distance = np.linalg.norm(origin - destination)

    assert_allclose(actual_distance, distances)
    assert (1, 1) == distances.shape


def test_calculate_points_within_distance(mock_case):

    origin = mock_case.electric.bipolar_egm.points
    destination = mock_case.points

    points_within_distance, distances = calculate_points_within_distance(
        origin,
        destination,
        max_distance=0,
        return_distances=True,
    )

    assert not (points_within_distance).all()
    assert 10 == np.sum(points_within_distance)

    assert_allclose(distances == 0, points_within_distance)


def test_calculate_points_within_distance_all(mock_case):

    origin = mock_case.electric.bipolar_egm.points
    destination = mock_case.points

    points_within_distance = calculate_points_within_distance(
        origin,
        destination,
        max_distance=np.inf,
        return_distances=False,
    )

    assert (points_within_distance).all()


@pytest.fixture(scope='module')
def interpolator(real_case):

    n_mapping_points = real_case.electric.bipolar_egm.points.shape[0]
    amplitudes = np.arange(n_mapping_points)

    return Interpolator(
        points=real_case.electric.bipolar_egm.points,
        field=amplitudes,
    )


def test_interpolator(real_case, interpolator):

    n_surface_points = real_case.points.shape[0]
    interpolated_field = interpolator(real_case.points)

    assert n_surface_points == interpolated_field.size


def test_interpolator_cutoff(real_case, interpolator):

    interpolated_field = interpolator(real_case.points, max_distance=0)

    assert interpolated_field.size == np.sum(np.isnan(interpolated_field))


def test_interpolator_kws(real_case):

    n_mapping_points = real_case.electric.bipolar_egm.points.shape[0]
    amplitudes = np.arange(n_mapping_points)

    kws_interpolator = Interpolator(
        points=real_case.electric.bipolar_egm.points,
        field=amplitudes,
        method_kws={'smoothing': 100},
    )

    assert 100 == kws_interpolator.method_kws['smoothing']


def test_interpolator_nearest(real_case):

    n_mapping_points = real_case.electric.bipolar_egm.points.shape[0]
    amplitudes = np.arange(n_mapping_points)

    nearest_interpolator = Interpolator(
        points=real_case.electric.bipolar_egm.points,
        field=amplitudes,
        method=scipy.interpolate.NearestNDInterpolator,
    )

    assert scipy.interpolate.NearestNDInterpolator is nearest_interpolator.method
    assert {} == nearest_interpolator.method_kws


def test_interpolate_voltage_onto_surface(real_case):

    # TODO: these regression tests could be refactored to use a mock data set

    n_surface_points = real_case.points.shape[0]
    interpolated_voltages = interpolate_voltage_onto_surface(real_case)

    assert n_surface_points == interpolated_voltages.size
    assert 5535 == np.sum(np.isnan(interpolated_voltages))


def test_interpolate_voltage_onto_surface_max_distance(real_case):

    n_surface_points = real_case.points.shape[0]
    interpolated_voltages = interpolate_voltage_onto_surface(real_case, max_distance=0)

    assert n_surface_points == np.sum(np.isnan(interpolated_voltages))
