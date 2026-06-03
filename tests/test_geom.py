#!/usr/bin/env python3
# coding: utf-8

import numpy as np
import pytest

from pdb_cpp import Coor, core, geom
from .datafiles import PDB_2OL9


def test_distance_matrix_ca_2ol9():
    coor = Coor(PDB_2OL9)
    ca = coor.select_atoms("name CA")

    dist = core.distance_matrix(ca.xyz, ca.xyz)

    print(dist)

    expected_dist = np.array(
        [
            [0.0, 3.802586, 7.114013, 10.46823, 13.52594, 16.591305],
            [3.802586, 0.0, 3.78049, 6.734785, 9.999435, 12.887392],
            [7.114013, 3.78049, 0.0, 3.8020015, 6.4243903, 9.637937],
            [10.46823, 6.734785, 3.8020015, 0.0, 3.8005443, 6.1563416],
            [13.525939, 9.999435, 6.4243903, 3.8005443, 0.0, 3.780037],
            [16.591307, 12.887392, 9.637937, 6.1563416, 3.780037, 0.0],
        ]
    )

    assert np.allclose(dist, expected_dist)

    assert dist.shape[0] == dist.shape[1] == ca.len
    assert np.allclose(np.diag(dist), 0.0)
    assert np.allclose(dist, dist.T)


def test_apply_transform_rowwise_coordinates():
    coords = np.array([[1.0, 2.0, 3.0], [0.0, 1.0, 0.0]])
    rotation = np.array(
        [
            [0.0, -1.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    translation = np.array([1.0, 2.0, 3.0])

    transformed = geom.apply_transform(coords, rotation, translation)
    expected = np.array([[-1.0, 3.0, 6.0], [0.0, 2.0, 3.0]])

    assert np.allclose(transformed, expected)


def test_xyz_setter_updates_coordinates_and_validates_shape():
    coor = Coor(PDB_2OL9)
    shifted = coor.xyz + np.array([1.0, -2.0, 0.5])

    coor.xyz = shifted

    assert np.allclose(coor.xyz, shifted)
    assert np.allclose(np.asarray(coor.x), shifted[:, 0])
    assert np.allclose(np.asarray(coor.y), shifted[:, 1])
    assert np.allclose(np.asarray(coor.z), shifted[:, 2])

    with pytest.raises(ValueError, match="shape"):
        coor.xyz = shifted[:, :2]
