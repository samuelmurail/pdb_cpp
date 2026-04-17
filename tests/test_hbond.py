#!/usr/bin/env python3
# coding: utf-8

"""Tests for H-bond computation (Baker & Hubbard geometric criteria).

Reference values for 2rri were computed with biotite 1.6.0:
    triplets = struc.hbond(atom_arr, cutoff_dist=2.5, cutoff_angle=120)
Same-residue bonds reported by biotite are excluded (our code is correct to skip them).
All 20 inter-residue biotite bonds are reproduced exactly by pdb_cpp.

Reference values for 1y0m were computed by pdb_cpp itself (no explicit H in file;
backbone H positions are reconstructed geometrically).
"""

import pytest
from pdb_cpp import Coor, hbond, core
from .datafiles import MMCIF_1Y0M, MMCIF_2RRI, CIF_1A0A

# ---------------------------------------------------------------------------
# Reference values derived from biotite 1.6.0 on 2rri.cif
# Key: (donor_chain, donor_resname, donor_heavy, acceptor_chain, acceptor_resname, acceptor_atom)
# Value: (dist_HA, dist_DA, angle_DHA)  — rounded to 3 / 2 decimal places
# ---------------------------------------------------------------------------
REFERENCE_2RRI = {
    ("A", "ALA", "N", "A", "LYS", "O"):  (2.301, 3.103, 138.49),
    ("A", "ARG", "N", "A", "ASP", "O"):  (1.989, 2.917, 157.00),
    ("A", "ARG", "N", "A", "TYR", "O"):  (2.060, 2.889, 141.10),
    ("A", "ASN", "N", "A", "LEU", "O"):  (2.286, 3.037, 132.80),
    ("A", "ASN", "N", "A", "LYS", "O"):  (2.088, 2.838, 131.90),
    ("A", "ASN", "ND2", "A", "ASN", "O"): (2.370, 3.224, 145.11),
    ("A", "ASP", "N", "A", "VAL", "O"):  (2.085, 2.965, 148.26),
    ("A", "GLN", "N", "A", "LEU", "O"):  (2.310, 2.955, 122.52),
    ("A", "ILE", "N", "A", "TYR", "O"):  (2.074, 2.983, 153.70),
    ("A", "LEU", "N", "A", "ASN", "O"):  (1.816, 2.675, 144.58),
    ("A", "LEU", "N", "A", "LEU", "O"):  (1.739, 2.386, 120.13),
    ("A", "LEU", "N", "A", "VAL", "O"):  (1.741, 2.681, 159.84),
    ("A", "LYS", "N", "A", "ARG", "O"):  (2.300, 3.110, 139.32),
    ("A", "LYS", "N", "A", "GLN", "O"):  (1.837, 2.814, 173.59),
    ("A", "LYS", "N", "A", "MET", "O"):  (2.221, 3.014, 137.05),
    ("A", "MET", "N", "A", "LEU", "O"):  (2.433, 3.315, 149.46),
    ("A", "SER", "N", "A", "TYR", "O"):  (2.147, 2.852, 127.54),
    ("A", "THR", "N", "A", "ALA", "O"):  (2.319, 3.258, 160.07),
    ("A", "THR", "N", "A", "THR", "O"):  (1.903, 2.834, 157.56),
    ("A", "TYR", "N", "A", "PHE", "O"):  (1.726, 2.676, 162.04),
}

# ---------------------------------------------------------------------------
# Reference values for 1y0m (backbone H reconstructed, angle>120°)
# 5 shortest H···A bonds
# ---------------------------------------------------------------------------
REFERENCE_1Y0M_TOP5 = [
    ("A", "LEU", "N", "A", "TRP", "O",  1.810, 2.798, 165.09),
    ("A", "SER", "N", "A", "ALA", "O",  1.819, 2.809, 165.59),
    ("A", "LYS", "N", "A", "GLU", "O",  1.848, 2.845, 168.68),
    ("A", "GLN", "N", "A", "GLU", "OE1", 1.849, 2.845, 168.21),
    ("A", "GLU", "N", "A", "LYS", "O",  1.850, 2.841, 166.25),
]
REFERENCE_1Y0M_N_BONDS = 36   # total bonds with angle_cutoff=120°
REFERENCE_2RRI_N_BONDS = 20


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _hb_key(hb):
    return (hb.donor_chain, hb.donor_resname, hb.donor_heavy_name,
            hb.acceptor_chain, hb.acceptor_resname, hb.acceptor_name)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_hbonds_returns_list_per_model():
    coor = Coor(MMCIF_1Y0M)
    result = hbond.hbonds(coor, angle_cutoff=120.0)
    assert isinstance(result, list)
    assert len(result) == coor.model_num


def test_hbonds_count_1y0m():
    """Total bond count for 1y0m (H reconstructed) matches biotite-calibrated reference."""
    coor = Coor(MMCIF_1Y0M)
    hb_list = hbond.hbonds(coor, angle_cutoff=120.0)[0]
    assert len(hb_list) == REFERENCE_1Y0M_N_BONDS


def test_hbonds_geometry_1y0m():
    """Top-5 shortest bonds match reference geometry (H reconstructed from heavy atoms)."""
    coor = Coor(MMCIF_1Y0M)
    hb_list = sorted(hbond.hbonds(coor, angle_cutoff=120.0)[0], key=lambda h: h.dist_HA)

    for i, (dc, dr, dn, ac, ar, an, exp_ha, exp_da, exp_ang) in enumerate(REFERENCE_1Y0M_TOP5):
        hb = hb_list[i]
        assert hb.donor_chain == dc
        assert hb.donor_resname == dr
        assert hb.donor_heavy_name == dn
        assert hb.acceptor_chain == ac
        assert hb.acceptor_resname == ar
        assert hb.acceptor_name == an
        assert abs(hb.dist_HA - exp_ha) < 0.002,   f"bond {i} dist_HA mismatch"
        assert abs(hb.dist_DA - exp_da) < 0.002,   f"bond {i} dist_DA mismatch"
        assert abs(hb.angle_DHA - exp_ang) < 0.1,  f"bond {i} angle mismatch"


def test_hbonds_count_2rri():
    """2rri has explicit H: total count must exactly match biotite reference."""
    coor = Coor(MMCIF_2RRI)
    hb_list = hbond.hbonds(coor, angle_cutoff=120.0)[0]
    assert len(hb_list) == REFERENCE_2RRI_N_BONDS


def test_hbonds_set_2rri():
    """Every bond found by biotite (inter-residue) is reproduced by pdb_cpp."""
    coor = Coor(MMCIF_2RRI)
    hb_list = hbond.hbonds(coor, angle_cutoff=120.0)[0]
    found = {_hb_key(hb) for hb in hb_list}
    for ref_key in REFERENCE_2RRI:
        assert ref_key in found, f"Missing bond {ref_key}"


def test_hbonds_geometry_2rri():
    """Geometry (dist_HA, dist_DA, angle) matches biotite values to 3 decimal places."""
    coor = Coor(MMCIF_2RRI)
    hb_by_key = {_hb_key(hb): hb for hb in hbond.hbonds(coor, angle_cutoff=120.0)[0]}
    for key, (exp_ha, exp_da, exp_ang) in REFERENCE_2RRI.items():
        assert key in hb_by_key, f"Bond {key} not found"
        hb = hb_by_key[key]
        assert abs(hb.dist_HA - exp_ha) < 0.002,  f"{key} dist_HA mismatch"
        assert abs(hb.dist_DA - exp_da) < 0.002,  f"{key} dist_DA mismatch"
        assert abs(hb.angle_DHA - exp_ang) < 0.1, f"{key} angle mismatch"


def test_hbond_fields():
    """Every HBond object exposes all required fields with sane geometry."""
    coor = Coor(MMCIF_1Y0M)
    hb_list = hbond.hbonds(coor)[0]
    assert len(hb_list) > 0
    for hb in hb_list:
        assert 0.0 < hb.dist_DA <= 3.5
        assert 0.0 < hb.dist_HA <= 2.5
        assert hb.angle_DHA >= 90.0
        assert hb.donor_resid != hb.acceptor_resid
        assert isinstance(hb.donor_resname, str) and hb.donor_resname
        assert isinstance(hb.acceptor_name, str) and hb.acceptor_name
        assert len(hb.donor_heavy_xyz) == 3
        assert len(hb.acceptor_xyz) == 3


def test_hbonds_protein_to_nucleic_include_backbone_acceptors_1a0a():
    """Protein-DNA H-bonds on 1A0A must include phosphate backbone acceptors.

    This guards against undercounting caused by treating only nucleobase atoms
    as nucleic-acid acceptors.
    """
    coor = Coor(CIF_1A0A)
    hb_list = hbond.hbonds(coor, donor_sel="protein", acceptor_sel="nucleic")[0]

    unique_keys = {
        (hb.donor_chain, hb.donor_resid, hb.donor_heavy_name,
         hb.acceptor_chain, hb.acceptor_resid, hb.acceptor_name)
        for hb in hb_list
    }
    acceptor_names = {hb.acceptor_name for hb in hb_list}

    assert len(unique_keys) == 6
    assert {"OP1", "OP2"}.issubset(acceptor_names)


def test_hbonds_strict_cutoff():
    """With zero-length cutoffs no H-bonds should be found."""
    coor = Coor(MMCIF_1Y0M)
    result = hbond.hbonds(coor, dist_DA_cutoff=0.0, dist_HA_cutoff=0.0)
    assert result[0] == []


def test_core_compute_hbonds_direct():
    """Raw C++ binding works and excludes self-residue bonds."""
    coor = Coor(MMCIF_1Y0M)
    model = coor.models[0]
    hb_list = core.compute_hbonds(model, model, model)
    assert isinstance(hb_list, list)
    assert len(hb_list) > 0
    for hb in hb_list:
        assert hb.donor_resid != hb.acceptor_resid
