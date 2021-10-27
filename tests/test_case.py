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
