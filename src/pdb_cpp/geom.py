#!/usr/bin/env python3
# coding: utf-8

from .core import distance_matrix


distance_matrix.__doc__ = """Compute a pairwise distance matrix between coordinates.

Parameters
----------
xyz_a : array-like
	Array of shape (N, 3) with coordinates.
xyz_b : array-like
	Array of shape (M, 3) with coordinates.

Returns
-------
numpy.ndarray
	Distance matrix of shape (N, M).
"""
