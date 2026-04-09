#!/usr/bin/env python3
# coding: utf-8

"""Tests for PQR and GRO format support.

Ported from pdb_numpy's test_gro.py and test_pqr.py, adapted to pdb_cpp API.
"""

import os
import pytest
import numpy as np

from pdb_cpp import Coor
from .datafiles import PDB_1Y0M, MMCIF_1Y0M, PQR_1Y0M, GRO_1Y0M, GRO_2RRI


# ---------------------------------------------------------------------------
# GRO format  (mirrors pdb_numpy/tests/test_gro.py)
# ---------------------------------------------------------------------------


def test_read_gro(tmp_path):
    """Read a multi-model GRO file (2rri, 20 NMR models)."""
    test = Coor(GRO_2RRI)
    assert test.len == 479
    assert test.model_num == 20

    assert test.resname_str[0] == "HIS"
    assert test.resid[0] == 1
    assert test.uniq_resid[0] == 0
    assert test.name_str[0] == "N"
    assert test.num[0] == 1
    assert test.x[0] == pytest.approx(-11.43, abs=1e-2)
    assert test.y[0] == pytest.approx(14.76, abs=1e-2)
    assert test.z[0] == pytest.approx(-14.63, abs=1e-2)
    np.testing.assert_allclose(
        test.xyz[0, :], np.array([-11.43, 14.76, -14.63]), atol=1e-2
    )


def test_read_pdb_write_gro(tmp_path):
    """PDB 1y0m → GRO → read back, check properties and coordinates."""
    test = Coor(PDB_1Y0M)
    test.write(os.path.join(tmp_path, "test.gro"))

    test_gro = Coor(os.path.join(tmp_path, "test.gro"))

    assert test.len == test_gro.len
    assert test.model_num == test_gro.model_num

    assert test_gro.resname_str[0] == "THR"
    assert test_gro.resid[0] == 791
    assert test_gro.uniq_resid[0] == 0
    assert test_gro.name_str[0] == "N"
    assert test_gro.num[0] == 1
    assert test_gro.x[0] == pytest.approx(-1.43, abs=1e-2)
    assert test_gro.y[0] == pytest.approx(9.76, abs=1e-2)
    assert test_gro.z[0] == pytest.approx(11.44, abs=1e-2)

    np.testing.assert_allclose(test.xyz, test_gro.xyz, atol=1e-2)


def test_read_mmcif_write_gro(tmp_path):
    """mmCIF 1y0m → GRO → read back, check properties and coordinates."""
    test = Coor(MMCIF_1Y0M)
    assert test.len == 648

    assert test.resname_str[0] == "THR"
    assert test.resid[0] == 791
    assert test.uniq_resid[0] == 0
    assert test.name_str[0] == "N"
    assert test.num[0] == 1
    assert test.x[0] == pytest.approx(-1.43, abs=1e-2)
    assert test.y[0] == pytest.approx(9.76, abs=1e-2)
    assert test.z[0] == pytest.approx(11.44, abs=1e-2)

    test.write(os.path.join(tmp_path, "test.gro"))
    test_gro = Coor(os.path.join(tmp_path, "test.gro"))

    assert test_gro.resname_str[0] == "THR"
    assert test_gro.resid[0] == 791
    assert test_gro.uniq_resid[0] == 0
    assert test_gro.name_str[0] == "N"
    assert test_gro.num[0] == 1
    assert test_gro.x[0] == pytest.approx(-1.43, abs=1e-2)
    assert test_gro.y[0] == pytest.approx(9.76, abs=1e-2)
    assert test_gro.z[0] == pytest.approx(11.44, abs=1e-2)


def test_read_gro_write_mmcif(tmp_path):
    """GRO 2rri → mmCIF → read back, verify multi-model and coordinates."""
    test_gro = Coor(GRO_2RRI)
    test_gro.write(os.path.join(tmp_path, "test_2RRI.cif"))
    test_cif = Coor(os.path.join(tmp_path, "test_2RRI.cif"))

    assert test_cif.len == 479
    assert test_cif.model_num == 20

    assert test_cif.resname_str[0] == "HIS"
    assert test_cif.resid[0] == 1
    assert test_cif.uniq_resid[0] == 0
    assert test_cif.name_str[0] == "N"
    assert test_cif.num[0] == 1
    assert test_cif.x[0] == pytest.approx(-11.43, abs=1e-2)
    assert test_cif.y[0] == pytest.approx(14.76, abs=1e-2)
    assert test_cif.z[0] == pytest.approx(-14.63, abs=1e-2)
    np.testing.assert_allclose(
        test_cif.xyz[0, :], np.array([-11.43, 14.76, -14.63]), atol=1e-2
    )

    np.testing.assert_allclose(test_cif.xyz, test_gro.xyz, atol=1e-2)


# ---------------------------------------------------------------------------
# PQR format  (mirrors pdb_numpy/tests/test_pqr.py)
# ---------------------------------------------------------------------------


def test_read_write_pqr(tmp_path):
    """Read PQR (with hydrogens), write PQR+PDB, compare non-H atoms."""
    test_pqr = Coor(PQR_1Y0M)
    assert test_pqr.len == 1362

    test_pqr.write(os.path.join(tmp_path, "test_2.pqr"))
    test_pqr.write(os.path.join(tmp_path, "test_2.pdb"))

    test_pdb = Coor(PDB_1Y0M)

    test_pdb_no_altloc = test_pdb.select_atoms("not altloc B C D")
    test_pqr_noh = test_pqr.select_atoms("not name H*")

    assert test_pqr_noh.len == test_pdb_no_altloc.len
    # Can't test the whole file because of the different atom order
    # in C-ter part
    np.testing.assert_allclose(
        test_pqr_noh.xyz[:500, :],
        test_pdb_no_altloc.xyz[:500, :],
        atol=1e-3,
    )


# ---------------------------------------------------------------------------
# Additional format tests
# ---------------------------------------------------------------------------


def test_gro_no_chain():
    """GRO has no chain IDs — chain defaults to '?'."""
    coor = Coor(GRO_1Y0M)
    assert coor.chain_str[0] == "?"


def test_pqr_selection():
    """Selections work on PQR-loaded structures."""
    coor = Coor(PQR_1Y0M)
    ca = coor.select_atoms("name CA")
    assert ca.len > 0
    assert ca.len < coor.len


def test_pdb_to_gro_to_pdb(tmp_path):
    """PDB → GRO → PDB roundtrip preserves structure."""
    pdb = Coor(PDB_1Y0M)
    gro_path = str(tmp_path / "step.gro")
    pdb_path = str(tmp_path / "step.pdb")
    pdb.write(gro_path)

    gro = Coor(gro_path)
    gro.write(pdb_path)

    pdb2 = Coor(pdb_path)
    assert pdb2.len == pdb.len
    np.testing.assert_allclose(pdb.xyz, pdb2.xyz, atol=0.01)


def test_pdb_to_pqr_to_pdb(tmp_path):
    """PDB → PQR → PDB roundtrip preserves coordinates."""
    pdb = Coor(PDB_1Y0M)
    pqr_path = str(tmp_path / "step.pqr")
    pdb_path = str(tmp_path / "step.pdb")
    pdb.write(pqr_path)

    pqr = Coor(pqr_path)
    pqr.write(pdb_path)

    pdb2 = Coor(pdb_path)
    assert pdb2.len == pdb.len
    np.testing.assert_allclose(pdb.xyz, pdb2.xyz, atol=0.01)


# ---------------------------------------------------------------------------
# format= parameter (explicit format override)
# ---------------------------------------------------------------------------


def test_read_format_pdb_explicit(tmp_path):
    """read() with format='pdb' parses a .pdb file correctly."""
    coor = Coor()
    coor.read(PDB_1Y0M, format="pdb")
    assert coor.len == 648
    assert coor.resname_str[0] == "THR"


def test_read_format_cif_explicit(tmp_path):
    """read() with format='cif' parses a .cif file correctly."""
    coor = Coor()
    coor.read(MMCIF_1Y0M, format="cif")
    assert coor.len == 648
    assert coor.resname_str[0] == "THR"


def test_read_format_pdb_wrong_extension(tmp_path):
    """format= overrides a misleading file extension (PDB content, .txt name)."""
    src = Coor(PDB_1Y0M)
    fake_path = str(tmp_path / "structure.txt")
    import shutil
    shutil.copy(PDB_1Y0M, fake_path)

    coor = Coor()
    assert coor.read(fake_path, format="pdb")
    assert coor.len == 648


def test_read_format_cif_wrong_extension(tmp_path):
    """format= overrides a misleading file extension (mmCIF content, .txt name)."""
    fake_path = str(tmp_path / "structure.txt")
    import shutil
    shutil.copy(MMCIF_1Y0M, fake_path)

    coor = Coor()
    assert coor.read(fake_path, format="cif")
    assert coor.len == 648


def test_read_format_pqr_explicit():
    """read() with format='pqr' parses a .pqr file correctly."""
    coor = Coor()
    coor.read(PQR_1Y0M, format="pqr")
    assert coor.len == 1362


def test_read_format_gro_explicit():
    """read() with format='gro' parses a .gro file correctly."""
    coor = Coor()
    coor.read(GRO_1Y0M, format="gro")
    assert coor.len == 648


def test_read_format_unknown_returns_false(tmp_path):
    """read() with an unknown format returns False without raising."""
    coor = Coor()
    result = coor.read(PDB_1Y0M, format="xyz")
    assert result is False
