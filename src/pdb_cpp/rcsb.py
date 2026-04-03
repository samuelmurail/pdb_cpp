#!/usr/bin/env python3
# coding: utf-8

"""Helpers for downloading and loading structures from the RCSB PDB."""

import os
import tempfile
import urllib.request
from urllib.error import HTTPError, URLError

from .core import Coor

__all__ = [
    "build_download_url",
    "download_structure",
    "load_structure",
    "download",
    "load",
]


_RCSB_DOWNLOAD_BASE = "https://files.rcsb.org/download"
_STRUCTURE_ALIASES = {
    "asymmetric_unit": "asymmetric_unit",
    "asym_unit": "asymmetric_unit",
    "asym": "asymmetric_unit",
    "entry": "asymmetric_unit",
    "model": "asymmetric_unit",
    "deposited": "asymmetric_unit",
    "biological_assembly": "biological_assembly",
    "biological assembly": "biological_assembly",
    "bioassembly": "biological_assembly",
    "assembly": "biological_assembly",
    "biounit": "biological_assembly",
}
_FORMAT_ALIASES = {
    "cif": "cif",
    "mmcif": "cif",
    "pdbx": "cif",
    "pdb": "pdb",
}


def _normalize_pdb_id(pdb_id):
    pdb_id = str(pdb_id).strip().lower()
    if not pdb_id:
        raise ValueError("pdb_id must be a non-empty string")
    return pdb_id


def _normalize_structure(structure):
    structure_key = str(structure).strip().lower().replace("-", "_")
    try:
        return _STRUCTURE_ALIASES[structure_key]
    except KeyError as exc:
        raise ValueError(
            "structure must be one of: asymmetric_unit, biological_assembly"
        ) from exc


def _normalize_file_format(file_format):
    format_key = str(file_format).strip().lower()
    try:
        return _FORMAT_ALIASES[format_key]
    except KeyError as exc:
        raise ValueError("file_format must be one of: cif, pdb") from exc


def _normalize_assembly_id(assembly_id):
    assembly_id = int(assembly_id)
    if assembly_id < 1:
        raise ValueError("assembly_id must be greater than or equal to 1")
    return assembly_id


def _get_cache_dir(cache_dir=None):
    if cache_dir is None:
        cache_dir = os.path.join(tempfile.gettempdir(), "pdb_cpp_cache", "rcsb")
    os.makedirs(cache_dir, mode=0o700, exist_ok=True)
    return cache_dir


def _get_local_filename(pdb_id, structure, file_format, assembly_id):
    if structure == "asymmetric_unit":
        return f"{pdb_id}.{file_format}"
    return f"{pdb_id}-assembly{assembly_id}.{file_format}"


def build_download_url(
    pdb_id,
    structure="asymmetric_unit",
    file_format="cif",
    assembly_id=1,
):
    """Build an RCSB download URL for a structure file.

    Parameters
    ----------
    pdb_id : str
        PDB identifier.
    structure : str, default="asymmetric_unit"
        Either the deposited asymmetric unit or a biological assembly.
    file_format : str, default="cif"
        Download format. Supported values are ``"cif"`` and ``"pdb"``.
    assembly_id : int, default=1
        Biological assembly identifier when ``structure`` is
        ``"biological_assembly"``.

    Returns
    -------
    str
        Download URL.
    """
    pdb_id = _normalize_pdb_id(pdb_id)
    structure = _normalize_structure(structure)
    file_format = _normalize_file_format(file_format)

    if structure == "asymmetric_unit":
        return f"{_RCSB_DOWNLOAD_BASE}/{pdb_id}.{file_format}"

    assembly_id = _normalize_assembly_id(assembly_id)
    if file_format == "cif":
        return f"{_RCSB_DOWNLOAD_BASE}/{pdb_id}-assembly{assembly_id}.cif"
    return f"{_RCSB_DOWNLOAD_BASE}/{pdb_id}.pdb{assembly_id}"


def download_structure(
    pdb_id,
    structure="asymmetric_unit",
    file_format="cif",
    assembly_id=1,
    cache_dir=None,
    force_download=False,
):
    """Download and cache an RCSB structure file.

    Parameters
    ----------
    pdb_id : str
        PDB identifier.
    structure : str, default="asymmetric_unit"
        Either ``"asymmetric_unit"`` or ``"biological_assembly"``.
    file_format : str, default="cif"
        Download format. Supported values are ``"cif"`` and ``"pdb"``.
    assembly_id : int, default=1
        Assembly identifier for biological assemblies.
    cache_dir : str, optional
        Cache directory. Defaults to a temporary directory managed by
        ``pdb_cpp``.
    force_download : bool, default=False
        Re-download the file even when it is already cached.

    Returns
    -------
    str
        Local path to the cached structure file.
    """
    pdb_id = _normalize_pdb_id(pdb_id)
    structure = _normalize_structure(structure)
    file_format = _normalize_file_format(file_format)
    assembly_id = _normalize_assembly_id(assembly_id)
    cache_dir = _get_cache_dir(cache_dir)

    local_name = _get_local_filename(pdb_id, structure, file_format, assembly_id)
    local_path = os.path.join(cache_dir, local_name)
    if os.path.exists(local_path) and not force_download:
        return local_path

    url = build_download_url(
        pdb_id,
        structure=structure,
        file_format=file_format,
        assembly_id=assembly_id,
    )
    try:
        with urllib.request.urlopen(url) as response:
            data = response.read()
    except (HTTPError, URLError) as exc:
        details = f"assembly {assembly_id}" if structure == "biological_assembly" else structure
        raise ValueError(
            f"Failed to fetch {details} for PDB ID {pdb_id} from {url}"
        ) from exc

    with open(local_path, "wb") as handle:
        handle.write(data)
    return local_path


def load_structure(
    pdb_id,
    structure="asymmetric_unit",
    file_format="cif",
    assembly_id=1,
    cache_dir=None,
    force_download=False,
):
    """Download a structure from RCSB and return it as a ``Coor`` object."""
    local_path = download_structure(
        pdb_id,
        structure=structure,
        file_format=file_format,
        assembly_id=assembly_id,
        cache_dir=cache_dir,
        force_download=force_download,
    )
    return Coor(local_path)


def download(*args, **kwargs):
    """Alias for :func:`download_structure`."""
    return download_structure(*args, **kwargs)


def load(*args, **kwargs):
    """Alias for :func:`load_structure`."""
    return load_structure(*args, **kwargs)