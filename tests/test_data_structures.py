import pytest
from numpy.testing import assert_allclose

import numpy as np
import pyvista

from openep.data_structures.case import Case
from openep.data_structures.surface import Fields
from openep._datasets.simple_meshes import CUBE

@pytest.fixture(scope='module')
def mesh():
    return pyvista.read(CUBE)

@pytest.fixture(scope='module')
def fields():
    
    data = np.arange(5)
    return Fields(
        data,
        data * 2,
        data * 3,
        data * 4,
        data * 5,
    )

@pytest.fixture(scope='module')
def case(mesh, fields):

    name = "Pretend-Case"

    points = mesh.points
    indices = mesh.faces.reshape(-1, 4)[:, 1:]  # ignore the number of vertices per face

    electric = None
    ablation = None
    notes = ['Note one', 'Note two']

    return Case(
        name=name,
        points=points,
        indices=indices,
        fields=fields,
        electric=electric,
        ablation=ablation,
        notes=notes,
    )

def test_case_creation(case):

    attributes = ['points', 'indices', 'fields', 'electric', 'ablation', 'notes']

    assert isinstance(case, Case)
    for attribute in attributes:
        assert hasattr(case, attribute)

def test_fields_access_with_key(case):
    
    data = np.arange(5)

    assert_allclose(data, case.fields.bipolar_voltage)

def test_mesh_creation(mesh, case):

    mesh_from_case = case.create_mesh()

    assert isinstance(mesh, pyvista.PolyData)
    assert mesh == mesh_from_case

def test_mesh_creation_back_faces(mesh, case):

    mesh_from_case = case.create_mesh(back_faces=True)
    assert mesh.n_faces == mesh_from_case.n_faces / 2

def test_mesh_creation_not_centered(mesh, case):

    mesh_from_case = case.create_mesh(recenter=False)

    assert_allclose(0, mesh.center)
    not assert_allclose(mesh.center, mesh_from_case.center)

def test_get_surface_data(case):

    points, indices = case.get_surface_data()
    assert case.points is points
    assert case.indices is indices

def test_get_surface_data_copy(case):

    points, indices = case.get_surface_data(copy=True)

    assert case.points is not points
    assert case.indices is not indices

    assert_allclose(points, case.points)
    assert_allclose(indices, case.indices)

def test_get_field(case):

    field = case.get_field(fieldname='force')
    assert case.fields['force'] is field

def test_get_field_copy(case):

    field = case.get_field(fieldname='force', copy=True)

    assert case.fields['force'] is not field
    assert_allclose(case.fields['force'], field)

def test_check_field(case):

    field = 'bipolar_voltage'
    missing_field = 'bubblegum'

    assert field in case.fields
    assert missing_field not in case.fields

def test_no_field(case):

    missing_field = 'bubblegum'
    match = f"There is no field '{missing_field}'"
    with pytest.raises(ValueError, match=match):
        case.fields[missing_field]

