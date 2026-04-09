#!/usr/bin/env python3
# coding: utf-8

"""DockQ reference-value regression tests."""

from .datafiles import (
    DOCKQ_MODEL,
    DOCKQ_NATIVE,
    PDB_1JD4,
    PDB_1RXZ,
    PDB_1RXZ_Colabfold,
    PDB_5M6N,
    PDB_1A2K,
    PDB_1A2K_MODEL,
    PDB_DIMER_DIMER,
    PDB_DIMER_DIMER_MODEL,
    CIF_1A0A,
    CIF_FOLD_2026_DNAPROT_MODEL,
)
from .dockq_reference_values import DOCKQ_REFERENCES, DOCKQ_MULTIMER_REFERENCES
from pdb_cpp import Coor, analysis
import pytest


def test_dockq_reference_1rxz_colabfold_vs_native():
    model_coor = Coor(PDB_1RXZ_Colabfold)
    native_coor = Coor(PDB_1RXZ)

    result = analysis.dockQ(model_coor, native_coor)
    ref = DOCKQ_REFERENCES["1rxz_colabfold_vs_1rxz"]

    assert result["DockQ"][0] == pytest.approx(ref["DockQ"], abs=0.001)
    assert result["Fnat"][0] == pytest.approx(ref["Fnat"], abs=0.001)
    assert result["Fnonnat"][0] == pytest.approx(ref["Fnonnat"], abs=0.001)
    assert result["iRMS"][0] == pytest.approx(ref["iRMS"], abs=0.001)
    assert result["LRMS"][0] == pytest.approx(ref["LRMS"], abs=0.001)


def test_dockq_reference_model_vs_native():
    model_coor = Coor(DOCKQ_MODEL)
    native_coor = Coor(DOCKQ_NATIVE)

    result = analysis.dockQ(model_coor, native_coor)
    ref = DOCKQ_REFERENCES["model_vs_native"]

    assert result["DockQ"][0] == pytest.approx(ref["DockQ"], abs=0.001)
    assert result["Fnat"][0] == pytest.approx(ref["Fnat"], abs=0.001)
    assert result["Fnonnat"][0] == pytest.approx(ref["Fnonnat"], abs=0.001)
    assert result["iRMS"][0] == pytest.approx(ref["iRMS"], abs=0.001)
    assert result["LRMS"][0] == pytest.approx(ref["LRMS"], abs=0.001)


def test_dockq_reference_1jd4_vs_5m6n_score_target():
    model_coor = Coor(PDB_1JD4)
    native_coor = Coor(PDB_5M6N)

    result = analysis.dockQ(model_coor, native_coor)
    ref = DOCKQ_REFERENCES["1jd4_vs_5m6n"]

    assert result["DockQ"][0] == pytest.approx(ref["DockQ"], abs=0.005)
    assert result["Fnat"][0] == pytest.approx(ref["Fnat"], abs=0.001)
    assert result["Fnonnat"][0] == pytest.approx(ref["Fnonnat"], abs=0.001)


def test_dockq_multimer_dimer_dimer_perfect():
    """4-chain perfect model (model == native): GlobalDockQ must be 1.000."""
    model_coor = Coor(PDB_DIMER_DIMER_MODEL)
    native_coor = Coor(PDB_DIMER_DIMER)
    ref = DOCKQ_MULTIMER_REFERENCES["dimer_dimer"]

    result = analysis.dockQ_multimer(model_coor, native_coor)

    # Correct number of interfaces with actual contacts
    valid_ifaces = [
        k for k, v in result["interfaces"].items()
        if v is not None and v["iRMS"][0] is not None
    ]
    assert len(valid_ifaces) == ref["n_interfaces"]

    # GlobalDockQ
    assert result["GlobalDockQ"][0] == pytest.approx(ref["GlobalDockQ"], abs=0.001)

    # Per-interface checks
    for iface_key, iref in ref["interfaces"].items():
        iface = result["interfaces"][iface_key]
        assert iface is not None
        assert iface["DockQ"][0] == pytest.approx(iref["DockQ"], abs=0.001)
        assert iface["Fnat"][0] == pytest.approx(iref["Fnat"], abs=0.001)


def test_dockq_multimer_1a2k_with_chain_map():
    """3-chain 1A2K with explicit optimal chain mapping (BAC:ABC)."""
    model_coor = Coor(PDB_1A2K_MODEL)
    native_coor = Coor(PDB_1A2K)
    # Optimal mapping found by external DockQ: native A→model B, B→A, C→C
    chain_map = {"A": "B", "B": "A", "C": "C"}
    ref = DOCKQ_MULTIMER_REFERENCES["1a2k_BAC_ABC"]

    result = analysis.dockQ_multimer(model_coor, native_coor, chain_map=chain_map)

    valid_ifaces = [
        k for k, v in result["interfaces"].items()
        if v is not None and v["iRMS"][0] is not None
    ]
    assert len(valid_ifaces) == ref["n_interfaces"]

    assert result["GlobalDockQ"][0] == pytest.approx(ref["GlobalDockQ"], abs=0.005)

    for iface_key, iref in ref["interfaces"].items():
        iface = result["interfaces"][iface_key]
        assert iface is not None, f"Interface {iface_key} missing from results"
        assert iface["DockQ"][0] == pytest.approx(iref["DockQ"], abs=0.005)
        assert iface["iRMS"][0] == pytest.approx(iref["iRMS"], abs=0.005)
        assert iface["LRMS"][0] == pytest.approx(iref["LRMS"], abs=0.005)
        assert iface["Fnat"][0] == pytest.approx(iref["Fnat"], abs=0.005)


def test_dockq_multimer_1a0a_protein_only():
    """Protein-DNA complex: only protein-protein interface scored.

    1A0A.cif has label chains A,B=DNA and C,D=protein (+water E-H).
    fold_2026 model has label chains A,B=protein and C,D=DNA.
    pdb_cpp does not yet support nucleic-acid interfaces; this test exercises
    dockQ_multimer with an explicit chain_map restricted to protein chains.
    External DockQ v2 GlobalDockQ over all 6 interfaces (protein+DNA) is 0.878.
    """
    model_coor = Coor(CIF_FOLD_2026_DNAPROT_MODEL)
    native_coor = Coor(CIF_1A0A)
    # native label C,D = protein; model label A,B = protein
    chain_map = {"C": "A", "D": "B"}
    ref = DOCKQ_MULTIMER_REFERENCES["1a0a_fold2026_protein_only"]

    result = analysis.dockQ_multimer(model_coor, native_coor, chain_map=chain_map)

    valid_ifaces = [
        k for k, v in result["interfaces"].items()
        if v is not None and v["iRMS"][0] is not None
    ]
    assert len(valid_ifaces) == ref["n_interfaces"]

    assert result["GlobalDockQ"][0] == pytest.approx(ref["GlobalDockQ"], abs=0.005)

    iface = result["interfaces"][("C", "D")]
    assert iface is not None
    iref = ref["interfaces"][("C", "D")]
    assert iface["DockQ"][0] == pytest.approx(iref["DockQ"], abs=0.005)
    assert iface["iRMS"][0] == pytest.approx(iref["iRMS"], abs=0.005)
    assert iface["LRMS"][0] == pytest.approx(iref["LRMS"], abs=0.005)
    assert iface["Fnat"][0] == pytest.approx(iref["Fnat"], abs=0.005)


def test_dockq_output_includes_clashes():
    """dockQ result dict must include a 'clashes' key with non-negative counts."""
    model_coor = Coor(DOCKQ_MODEL)
    native_coor = Coor(DOCKQ_NATIVE)

    result = analysis.dockQ(model_coor, native_coor)

    assert "clashes" in result
    assert len(result["clashes"]) == len(result["DockQ"])
    assert all(c >= 0 for c in result["clashes"])


def test_dockq_multimer_1a2k_auto_mapping():
    """dockQ_multimer without chain_map must find the optimal A→B, B→A, C→C mapping."""
    model_coor = Coor(PDB_1A2K_MODEL)
    native_coor = Coor(PDB_1A2K)
    ref = DOCKQ_MULTIMER_REFERENCES["1a2k_auto_mapping"]

    result = analysis.dockQ_multimer(model_coor, native_coor)  # no chain_map

    # Returned chain_map must be the optimal permutation.
    assert result["chain_map"] == ref["chain_map"]

    valid_ifaces = [
        k for k, v in result["interfaces"].items()
        if v is not None and v["iRMS"][0] is not None
    ]
    assert len(valid_ifaces) == ref["n_interfaces"]

    assert result["GlobalDockQ"][0] == pytest.approx(ref["GlobalDockQ"], abs=0.005)

    for iface_key, iref in ref["interfaces"].items():
        iface = result["interfaces"][iface_key]
        assert iface is not None, f"Interface {iface_key} missing from results"
        assert iface["DockQ"][0] == pytest.approx(iref["DockQ"], abs=0.005)
        assert iface["iRMS"][0] == pytest.approx(iref["iRMS"], abs=0.005)
        assert iface["LRMS"][0] == pytest.approx(iref["LRMS"], abs=0.005)
        assert iface["Fnat"][0] == pytest.approx(iref["Fnat"], abs=0.005)


def test_dockq_multimer_returns_chain_map():
    """dockQ_multimer result must always include a 'chain_map' key."""
    model_coor = Coor(PDB_DIMER_DIMER_MODEL)
    native_coor = Coor(PDB_DIMER_DIMER)

    # With explicit chain_map
    explicit_map = {"A": "A", "B": "B", "H": "H", "L": "L"}
    result_explicit = analysis.dockQ_multimer(model_coor, native_coor, chain_map=explicit_map)
    assert result_explicit["chain_map"] == explicit_map

    # Without chain_map (auto-detected)
    result_auto = analysis.dockQ_multimer(model_coor, native_coor)
    assert "chain_map" in result_auto
    assert set(result_auto["chain_map"].keys()) == set(explicit_map.keys())