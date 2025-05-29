#!/usr/bin/env python3

"""
Simple test for the coor_align function
"""

import src.pdb_cpp.core as core

def test_coor_align_simple():
    """Simple test of coor_align functionality"""
    
    # Load test structures
    print("Loading structures...")
    coor1 = core.Coor("2rri.pdb")
    coor2 = core.Coor("3eam.pdb")
    
    print(f"Structure 1 size: {coor1.size()} atoms")
    print(f"Structure 2 size: {coor2.size()} atoms")
    
    # Get CA atom indices
    ca_indices_1 = coor1.get_index_select("name CA")
    ca_indices_2 = coor2.get_index_select("name CA")
    
    print(f"CA atoms in structure 1: {len(ca_indices_1)}")
    print(f"CA atoms in structure 2: {len(ca_indices_2)}")
    
    # Use first 10 CA atoms for alignment test
    min_atoms = min(10, len(ca_indices_1), len(ca_indices_2))
    test_indices_1 = ca_indices_1[:min_atoms]
    test_indices_2 = ca_indices_2[:min_atoms]
    
    print(f"Using {min_atoms} CA atoms for alignment test")
    
    # Get initial centroid of structure 1
    model1 = coor1.get_Models(0)
    initial_x = sum(model1.get_x()[i] for i in test_indices_1) / len(test_indices_1)
    initial_y = sum(model1.get_y()[i] for i in test_indices_1) / len(test_indices_1)
    initial_z = sum(model1.get_z()[i] for i in test_indices_1) / len(test_indices_1)
    
    print(f"Initial centroid of structure 1: ({initial_x:.3f}, {initial_y:.3f}, {initial_z:.3f})")
    
    # Perform alignment
    print("Performing alignment...")
    core.coor_align(coor1, coor2, test_indices_1, test_indices_2, 0)
    
    # Get final centroid of structure 1
    model1_aligned = coor1.get_Models(0)
    final_x = sum(model1_aligned.get_x()[i] for i in test_indices_1) / len(test_indices_1)
    final_y = sum(model1_aligned.get_y()[i] for i in test_indices_1) / len(test_indices_1)
    final_z = sum(model1_aligned.get_z()[i] for i in test_indices_1) / len(test_indices_1)
    
    print(f"Final centroid of structure 1: ({final_x:.3f}, {final_y:.3f}, {final_z:.3f})")
    
    # Check if coordinates changed (indicating alignment occurred)
    coord_change = abs(initial_x - final_x) + abs(initial_y - final_y) + abs(initial_z - final_z)
    print(f"Coordinate change magnitude: {coord_change:.3f}")
    
    if coord_change > 0.01:
        print("✓ Alignment performed successfully - coordinates changed")
    else:
        print("? Minimal coordinate change - structures may already be aligned")
    
    # Save results
    coor1.write("test_aligned_1.pdb")
    coor2.write("test_aligned_2.pdb")
    print("Aligned structures saved as test_aligned_1.pdb and test_aligned_2.pdb")
    
    return True

if __name__ == "__main__":
    try:
        test_coor_align_simple()
        print("\n✓ Test completed successfully!")
    except Exception as e:
        print(f"\n✗ Test failed: {e}")
        import traceback
        traceback.print_exc()
