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
from numpy.testing import assert_allclose

import numpy as np

import openep
from openep._datasets.openep_datasets import DATASET_2


@pytest.fixture()
def case():
    return openep.load_openep_mat(DATASET_2)


@pytest.fixture()
def exported_case(case, tmp_path):

    EXPORTED_CASE = tmp_path / "exported_case.mat"
    openep.export_openep_mat(case, EXPORTED_CASE.as_posix())

    return openep.load_openep_mat(EXPORTED_CASE)


def test_openep_mat_export(case, exported_case):
    """Check the scalar fields of the original and exported data set are equal."""

    assert_allclose(case.fields.bipolar_voltage, exported_case.fields.bipolar_voltage, equal_nan=True)
    assert_allclose(case.fields.unipolar_voltage, exported_case.fields.unipolar_voltage, equal_nan=True)
    assert_allclose(case.fields.local_activation_time, exported_case.fields.local_activation_time, equal_nan=True)
    assert_allclose(case.fields.force, exported_case.fields.force, equal_nan=True)
    assert_allclose(case.fields.impedance, exported_case.fields.impedance, equal_nan=True)
