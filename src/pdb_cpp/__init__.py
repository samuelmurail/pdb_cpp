#!/usr/bin/env python3
# coding: utf-8

# Autorship information
__author__ = "Samuel Murail"
__copyright__ = "Copyright 2025, RPBS"
__credits__ = ["Samuel Murail"]
__license__ = "GNU General Public License v2.0"
__version__ = "0.0.1"
__maintainer__ = "Samuel Murail"
__email__ = "samuel.murail@u-paris.fr"
__status__ = "Beta"

from .core import Coor, Model

# Importing _pyprops applies runtime patches for python-facing helpers.
from . import _pyprops  # noqa: F401

__all__ = ["Coor", "Model"]
