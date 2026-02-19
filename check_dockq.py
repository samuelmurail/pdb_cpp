#!/usr/bin/env python3
# coding: utf-8

"""Run pdb_cpp dockQ with explicit chains to match DockQ CLI mapping."""

from pdb_cpp import Coor, analysis

# Replace with your file paths
model_path = "/home/murail/owncloud/Private/test/7AAM_unrelaxed_rank_049_alphafold2_multimer_v3_model_2_seed_012.pdb"
native_path = "/home/murail/owncloud/Private/test/7AAM.pdb"

model = Coor(model_path)
native = Coor(native_path)

# Explicit chain mapping: model A,B -> native A,C (AB:AC)
result = analysis.dockQ(
    model,
    native,
    rec_chains=['A'],      # Model receptor chain
    lig_chains=['B'],      # Model ligand chain
    native_rec_chains=['A'],  # Native receptor chain
    native_lig_chains=['C']   # Native ligand chain
)

print("pdb_cpp DockQ results with AB:AC mapping:")
print(f"DockQ: {result['DockQ'][0]:.3f}")
print(f"Fnat: {result['Fnat'][0]:.3f}")
print(f"Fnonnat: {result['Fnonnat'][0]:.3f}")
print(f"iRMS: {result['iRMS'][0]:.3f}")
print(f"LRMS: {result['LRMS'][0]:.3f}")

# Compare with automatic chain selection
result_auto = analysis.dockQ(model, native)
print("\nAutomatic chain selection:")
print(f"DockQ: {result_auto['DockQ'][0]:.3f}")
print(f"Fnat: {result_auto['Fnat'][0]:.3f}")
print(f"Fnonnat: {result_auto['Fnonnat'][0]:.3f}")
print(f"iRMS: {result_auto['iRMS'][0]:.3f}")
print(f"LRMS: {result_auto['LRMS'][0]:.3f}")