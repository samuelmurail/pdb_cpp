#!/usr/bin/env python3
# coding: utf-8

import numpy as np

from .core import distance_matrix, compute_dihedrals


def apply_transform(coords, rotation, translation):
	"""Apply a rigid transform to 3D coordinates.

	The transform follows $x' = R x + t$ with row-wise coordinate arrays,
	implemented as ``coords @ R.T + t``.

	Parameters
	----------
	coords : array-like, shape (N, 3)
		Input coordinates.
	rotation : array-like, shape (3, 3)
		Rotation matrix.
	translation : array-like, shape (3,)
		Translation vector.

	Returns
	-------
	numpy.ndarray
		Transformed coordinates with shape ``(N, 3)``.
	"""
	xyz = np.asarray(coords, dtype=float)
	if xyz.ndim != 2 or xyz.shape[1] != 3:
		raise ValueError("coords must have shape (N, 3)")

	rot = np.asarray(rotation, dtype=float)
	if rot.shape != (3, 3):
		raise ValueError("rotation must have shape (3, 3)")

	trans = np.asarray(translation, dtype=float)
	if trans.shape != (3,):
		raise ValueError("translation must have shape (3,)")

	return xyz @ rot.T + trans


__all__ = ["distance_matrix", "compute_dihedrals", "apply_transform"]
