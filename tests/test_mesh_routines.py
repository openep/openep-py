import pytest
from numpy.testing import assert_allclose

import numpy as np
import pyvista

from openep.mesh.mesh_routines import (
    calculate_field_area, calculate_mesh_volume,
    calculate_vertex_distance, calculate_vertex_path
)
from openep._datasets.simple_meshes import (
    CUBE, SPHERE, BROKEN_SPHERE, TRIANGLES
)


@pytest.fixture(scope='module')
def cube():
    return pyvista.read(CUBE)


@pytest.fixture(scope='module')
def sphere():
    return pyvista.read(SPHERE)


@pytest.fixture(scope='module')
def sphere_data(sphere):

    faces = sphere.faces.reshape(-1, 4)[:, 1:]
    triangles = sphere.points[faces]
    areas = sphere.compute_cell_sizes(
            length=False,
            area=True,
            volume=False,
        )['Area']

    return {
        "faces": faces,
        "triangles": triangles,
        "areas": areas,
    }


@pytest.fixture(scope='module')
def broken_sphere():
    return pyvista.read(BROKEN_SPHERE)


@pytest.fixture(scope='module')
def triangles():
    return pyvista.read(TRIANGLES)


def test_calculate_mesh_volume(sphere):

    volume = calculate_mesh_volume(sphere, fill_holes=False)
    assert_allclose(sphere.volume, volume)


def test_calculate_mesh_volume_repair(broken_sphere, sphere):

    volume = calculate_mesh_volume(broken_sphere, fill_holes=True)
    assert_allclose(sphere.volume, volume, atol=0.002)


def test_calculate_field_area(sphere, sphere_data):

    xbelow0 = sphere_data['triangles'][..., 0].mean(axis=1) <= 0  # select every triangle with mean x below YZ plane
    calculated_area = sphere_data['areas'][xbelow0].sum()

    area = calculate_field_area(sphere, sphere.points[:, 1], threshold=0)

    assert_allclose(calculated_area, area)


def test_calculate_vertex_distance_euclidian(cube):

    test_dist = calculate_vertex_distance(
        mesh=cube,
        start_index=0,
        end_index=7,
        metric='euclidian'
    )  # far corners are at indices 0 and 7

    assert_allclose(test_dist, np.sqrt(3), atol=5)


def test_calculate_vertex_distance_euclidian_sphere(sphere):

    sphere_diameter = 2.0
    start_index = 18  # top of sphere
    end_index = 23  # bottom of sphere
    test_dist = calculate_vertex_distance(
        mesh=sphere,
        start_index=start_index,
        end_index=end_index,
        metric='euclidian'
    )

    assert_allclose(test_dist, sphere_diameter, atol=5)


def test_calculate_vertex_distance_geodesic_sphere(sphere):

    sphere_half_circumference = 3.138363827815294  # should be equal to pi*diameter/2=pi, but it's not a perfect sphere
    start_index = 18  # top of sphere
    end_index = 23  # bottom of sphere
    test_dist = calculate_vertex_distance(
        mesh=sphere,
        start_index=start_index,
        end_index=end_index,
        metric='geodesic'
    )

    assert_allclose(test_dist, sphere_half_circumference, atol=5)


def test_calculate_vertex_path(cube):

    path = calculate_vertex_path(cube, 0, 7)
    assert_allclose(path, [0, 1, 7])
