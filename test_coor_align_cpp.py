#!/usr/bin/env python3

"""
Test the C++ coor_align function implementation
"""

import numpy as np
from src.pdb_cpp import Coor, core

def test_coor_align_basic():
    """Test basic coor_align functionality"""
    
    # Load test structures
    coor1 = Coor("5bkg.pdb")
    coor2 = Coor("3eam.pdb") 


    print(f"Coor1 size: {coor1.size()}")
    print(f"Coor2 size: {coor2.size()}")
    
    # Get CA atom indices for alignment
    ca_indices_1 = coor1.get_index_select("name CA")
    ca_indices_2 = coor2.get_index_select("name CA")
    
    print(f"CA atoms in coor1: {len(ca_indices_1)}")
    print(f"CA atoms in coor2: {len(ca_indices_2)}")
    
    common_atoms = core.get_common_atoms(coor1, coor2)
    if not common_atoms:
        raise ValueError("No common CA atoms found between the two structures.")
    print(f"Common CA atoms: {len(common_atoms[0])}")
    ca_indices_1 = common_atoms[0]
    ca_indices_2 = common_atoms[1]

    # For testing, use first min(len1, len2) atoms
    min_atoms = min(len(ca_indices_1), len(ca_indices_2))
    test_indices_1 = ca_indices_1[:min_atoms]
    test_indices_2 = ca_indices_2[:min_atoms]
    
    print(f"Using {min_atoms} CA atoms for alignment")
    
    # Get initial positions before alignment
    model1 = coor1.get_Models(0)
    model2 = coor2.get_Models(0)
    
    initial_pos_1 = [(model1.get_x()[i], model1.get_y()[i], model1.get_z()[i]) 
                     for i in test_indices_1[:5]]  # First 5 atoms
    initial_pos_2 = [(model2.get_x()[i], model2.get_y()[i], model2.get_z()[i]) 
                     for i in test_indices_2[:5]]  # First 5 atoms
    
    print("Initial positions (first 5 CA atoms):")
    print("Coor1:", initial_pos_1)
    print("Coor2:", initial_pos_2)
    
    # Perform alignment using C++ coor_align function
    print("\nPerforming C++ coor_align...")
    core.coor_align(coor1, coor2, test_indices_1, test_indices_2, 0)
    
    # Get positions after alignment
    model1_aligned = coor1.get_Models(0)
    model2_aligned = coor2.get_Models(0)
    
    aligned_pos_1 = [(model1_aligned.get_x()[i], model1_aligned.get_y()[i], model1_aligned.get_z()[i]) 
                     for i in test_indices_1[:5]]
    aligned_pos_2 = [(model2_aligned.get_x()[i], model2_aligned.get_y()[i], model2_aligned.get_z()[i]) 
                     for i in test_indices_2[:5]]
    
    print("Aligned positions (first 5 CA atoms):")
    print("Coor1:", aligned_pos_1)
    print("Coor2:", aligned_pos_2)
    
    # Calculate RMSD between aligned structures
    rmsd = 0.0
    for i in range(min_atoms):
        idx1, idx2 = test_indices_1[i], test_indices_2[i]
        dx = model1_aligned.get_x()[idx1] - model2_aligned.get_x()[idx2]
        dy = model1_aligned.get_y()[idx1] - model2_aligned.get_y()[idx2]
        dz = model1_aligned.get_z()[idx1] - model2_aligned.get_z()[idx2]
        rmsd += dx*dx + dy*dy + dz*dz
    
    rmsd = np.sqrt(rmsd / min_atoms)
    print(f"\nRMSD after alignment: {rmsd:.4f} Å")
    
    # Save aligned structures for verification
    coor1.write("tmp_cpp_aligned_1.pdb")
    coor2.write("tmp_cpp_aligned_2.pdb")
    print("Aligned structures saved as tmp_cpp_aligned_1.pdb and tmp_cpp_aligned_2.pdb")
    
    return rmsd

if __name__ == "__main__":
    try:
        rmsd = test_coor_align_basic()
        print(f"\n✓ Test completed successfully! Final RMSD: {rmsd:.4f} Å")
    except Exception as e:
        print(f"✗ Test failed with error: {e}")
        import traceback
        traceback.print_exc()
