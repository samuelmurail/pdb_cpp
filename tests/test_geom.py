#!/usr/bin/env python3
# coding: utf-8

import numpy as np

from pdb_cpp import Coor, core
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
