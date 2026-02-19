#!/usr/bin/env python3
# coding: utf-8

"""Run pdb_cpp dockQ with explicit chains to match DockQ CLI mapping."""

import logging
from pdb_cpp import Coor, analysis, core

logging.basicConfig(level=logging.INFO, format="%(levelname)s:%(name)s:%(message)s", force=True)
logging.getLogger("pdb_cpp").setLevel(logging.INFO)   # optional: only pdb_cpp logs


# Replace with your file paths
model_path = "/home/murail/owncloud/Private/test/7M6T_unrelaxed_rank_075_alphafold2_multimer_v3_model_1_seed_005.pdb"
native_path = "/home/murail/owncloud/Private/test/7M6T.pdb"

model = Coor(model_path)
native = Coor(native_path)

rec_chain_model = "A"
lig_chain_model = "B"
rec_chain_native = "A"
lig_chain_native = "D"


def common_backbone_atoms(model_chain: str, native_chain: str) -> int:
    model_idx, native_idx = core.get_common_atoms(
        model,
        native,
        chain_1=[model_chain],
        chain_2=[native_chain],
        back_names=["CA", "N", "C", "O"],
        matrix_file="",
    )
    return min(len(model_idx), len(native_idx))


native_chain_candidates = list(native.get_aa_seq().keys())
compat = {chain: common_backbone_atoms(lig_chain_model, chain) for chain in native_chain_candidates}

print("Ligand compatibility (model chain B vs native chains):")
for chain in native_chain_candidates:
    print(f"  B:{chain} -> {compat[chain]} common backbone atoms")

if compat.get(lig_chain_native, 0) == 0:
    suggested = max(compat, key=compat.get)
    print(
        f"WARNING: model chain {lig_chain_model} to native chain {lig_chain_native} has 0 common backbone atoms. "
        f"Using native ligand chain {suggested} instead ({compat[suggested]} common atoms)."
    )
    lig_chain_native = suggested

# Explicit chain mapping: model A,B -> native A,C (AB:AC)
result = analysis.dockQ(
    model,
    native,
    rec_chains=[rec_chain_model],
    lig_chains=[lig_chain_model],
    native_rec_chains=[rec_chain_native],
    native_lig_chains=[lig_chain_native],
)


def fmt_metric(value):
    if value is None:
        return "None"
    return f"{value:.3f}"

print("pdb_cpp DockQ results with AB:AC mapping:")
print(f"DockQ: {fmt_metric(result['DockQ'][0])}")
print(f"Fnat: {fmt_metric(result['Fnat'][0])}")
print(f"Fnonnat: {fmt_metric(result['Fnonnat'][0])}")
print(f"iRMS: {fmt_metric(result['iRMS'][0])}")
print(f"LRMS: {fmt_metric(result['LRMS'][0])}")

# Compare with automatic chain selection
print("\nAutomatic chain selection:")
try:
    result_auto = analysis.dockQ(model, native)
    print(f"DockQ: {fmt_metric(result_auto['DockQ'][0])}")
    print(f"Fnat: {fmt_metric(result_auto['Fnat'][0])}")
    print(f"Fnonnat: {fmt_metric(result_auto['Fnonnat'][0])}")
    print(f"iRMS: {fmt_metric(result_auto['iRMS'][0])}")
    print(f"LRMS: {fmt_metric(result_auto['LRMS'][0])}")
except ValueError as exc:
    print(f"Automatic mapping failed: {exc}")