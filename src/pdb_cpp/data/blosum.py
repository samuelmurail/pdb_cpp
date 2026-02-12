#!/usr/bin/env python3
# coding: utf-8

import os

_BLOSUM62_CACHE = None


def load_blosum(path=None):
    """Load a BLOSUM substitution matrix from a file.

    Parameters
    ----------
    path : str, optional
        Path to the matrix file. Defaults to blosum62.txt in this package.

    Returns
    -------
    dict
        Mapping of (aa1, aa2) to substitution score.
    """
    if path is None:
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "blosum62.txt")

    blosum = {}
    aa_list = []
    with open(path, "r") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            if line.startswith(" "):
                line = line.split()
                aa_list = line[:-1]
                continue

            line = line.split()
            for i, aa in enumerate(aa_list):
                blosum[(line[0], aa)] = int(line[i + 1])
                blosum[(aa, line[0])] = int(line[i + 1])

    return blosum


def get_blosum62():
    """Return a cached BLOSUM62 matrix.

    Returns
    -------
    dict
        Mapping of (aa1, aa2) to substitution score.
    """
    global _BLOSUM62_CACHE
    if _BLOSUM62_CACHE is None:
        _BLOSUM62_CACHE = load_blosum()
    return _BLOSUM62_CACHE


class _LazyBlosum62(dict):
    def _ensure(self):
        if not self:
            self.update(get_blosum62())

    def __getitem__(self, key):
        self._ensure()
        return super().__getitem__(key)

    def get(self, key, default=None):
        self._ensure()
        return super().get(key, default)

    def __contains__(self, key):
        self._ensure()
        return super().__contains__(key)


BLOSUM62 = _LazyBlosum62()

__all__ = ["BLOSUM62", "get_blosum62", "load_blosum"]
