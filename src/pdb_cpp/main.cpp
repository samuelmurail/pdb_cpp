#include <chrono>
#include <iomanip>

#include "Coor.h"
#include "select.h"
#include "Model.h"
#include "TMAlign.h"
#include "sequence_align.h"

using namespace std;

int main() {
    Coor structure;

    auto start = chrono::high_resolution_clock::now();
    structure.read("3eam_gap.pdb");
    //structure.read("2rri.pdb");
    auto end = chrono::high_resolution_clock::now();

    chrono::duration<double> elapsed = end - start;
    cout << "Time taken to get coordinates:   " << setprecision(3) << elapsed.count() << " seconds\n";

    // cout << "First atom: ";
    // for(int i = 0 ; i < 5 ; i ++ ){
    //     cout << structure.getAtomNames()[0][i] ;
    // }
    // cout << " ";
    // for(int i = 0 ; i < 5 ; i ++ ){
    //     cout << structure.getResNames()[0][i] ;
    // }
    // cout << " ";
    // // << structure.getAtomNames()[0] << " "
    // // << structure.getResNames()[0] << " "
    // cout << structure.getX()[0] << " "
    // << structure.getY()[0] << " "
    // << structure.getZ()[0] << endl;

    //structure.clear();


    start = chrono::high_resolution_clock::now();
    structure.write("tmp.pdb");
    end = chrono::high_resolution_clock::now();
    elapsed = end - start;
    cout << "Time taken to write coordinates: " << setprecision(3) << elapsed.count() << " seconds\n";

    //structure.transformation.print();
    //structure.symmetry.print();

    Model model = structure.get_Models(0);
    cout << "Model size: " << model.size() << endl;
    vector<bool> indexes = model.simple_select_atoms("name", {"CA"}, "isin");
    cout << "Number of CA atoms: " << count(indexes.begin(), indexes.end(), true) << endl;
    indexes = model.simple_select_atoms("name", {"CB"}, "isin");
    cout << "Number of CB atoms: " << count(indexes.begin(), indexes.end(), true) << endl;
    indexes = model.simple_select_atoms("name", {"CG"}, "==");
    cout << "Number of CG atoms: " << count(indexes.begin(), indexes.end(), true) << endl;
    start = chrono::high_resolution_clock::now();

    indexes = model.simple_select_atoms("name", {"CA", "CB", "CG"}, "isin");
    cout << "Number of CA CB CG atoms: " << count(indexes.begin(), indexes.end(), true) << endl;
    end = chrono::high_resolution_clock::now();
    elapsed = end - start;
    cout << "Time taken to select coordinates: " << setprecision(3) << elapsed.count() << " seconds\n";

    indexes = model.simple_select_atoms("resname", {"ALA"}, "==");
    cout << "Number of ALA atoms: " << count(indexes.begin(), indexes.end(), true) << endl;

    indexes = model.simple_select_atoms("resname", {"ALA", "GLY", "CYS"}, "isin");
    cout << "Number of ALA GLY CYS atoms: " << count(indexes.begin(), indexes.end(), true) << endl;
    end = chrono::high_resolution_clock::now();
    elapsed = end - start;
    cout << "Time taken to select coordinates: " << setprecision(3) << elapsed.count() << " seconds\n";

    indexes = model.simple_select_atoms("x", {"30.0"}, ">=");
    cout << "Number of x >= 0.0 atoms: " << count(indexes.begin(), indexes.end(), true) << endl;
    end = chrono::high_resolution_clock::now();
    elapsed = end - start;
    cout << "Time taken to select coordinates: " << setprecision(3) << elapsed.count() << " seconds\n";

    indexes = model.simple_select_atoms("z", {"20.0"}, "<");
    cout << "Number of z < 20.0 atoms: " << count(indexes.begin(), indexes.end(), true) << endl;

    indexes = model.simple_select_atoms("resid", {"20"}, "<");
    cout << "Number of resid < 20.0 atoms: " << count(indexes.begin(), indexes.end(), true) << endl;

    indexes = model.simple_select_atoms("resid", {"20", "21", "25"}, "isin");
    cout << "Number of resid 20 21 25 atoms: " << count(indexes.begin(), indexes.end(), true) << endl;

    indexes = model.simple_select_atoms("chain", {"A", "B"}, "isin");
    cout << "Number of chain A B atoms: " << count(indexes.begin(), indexes.end(), true) << endl;

    indexes = model.simple_select_atoms("name", {"C"}, "startswith");
    cout << "Number of name C* atoms: " << count(indexes.begin(), indexes.end(), true) << endl;

    string selection = "resname ALA GLY CYS and chain A B";
    Token parsed_selection = parse_selection(selection);
    cout << "Parsed selection: ";
    print_tokens(parsed_selection);
    cout << "Starting selection..." << endl;
    indexes = model.select_tokens(parsed_selection);
    cout << "Number of SEL atoms: " << count(indexes.begin(), indexes.end(), true) << endl;

    selection = "resid > 250";
    parsed_selection = parse_selection(selection);
    cout << "Parsed selection: ";
    print_tokens(parsed_selection);
    cout << "Starting selection..." << endl;
    indexes = model.select_tokens(parsed_selection);
    cout << "Number of SEL atoms: " << count(indexes.begin(), indexes.end(), true) << endl;

    selection = "resid > 250 and chain A B and resname ALA GLY";
    parsed_selection = parse_selection(selection);
    cout << "Parsed selection: ";
    print_tokens(parsed_selection);
    cout << "Starting selection..." << endl;
    indexes = model.select_tokens(parsed_selection);
    cout << "Number of SEL atoms: " << count(indexes.begin(), indexes.end(), true) << endl;

    indexes = model.select_atoms(selection);
    cout << "Number of SEL2 atoms: " << count(indexes.begin(), indexes.end(), true) << endl;

    selection = "resid >= 250 and not chain C D E and resname ALA GLY and within 20.0 of chain C";

    //parsed_selection = parse_selection(selection);
    //cout << "Parsed selection: ";
    //print_tokens(parsed_selection);
    cout << "Starting selection..." << endl;
    indexes = model.select_atoms(selection);
    cout << "- Number of not SEL atoms: " << count(indexes.begin(), indexes.end(), true) << endl;
    selection = "resid >= 250 and chain A B and resname ALA GLY and within 20.0 of chain C";
    indexes = model.select_atoms(selection);
    cout << "- Number of not SEL atoms: " << count(indexes.begin(), indexes.end(), true) << endl;

    vector<int> indices = model.get_index_select(selection);
    cout << "Number of indices: " << indices.size() << endl;
    cout << "First index: " << indices[0] << endl;
    cout << "Last index: " << indices[indices.size() - 1] << endl;
    cout << "Indexes: ";
    for (const auto &index : indices) {
        cout << index << " ";
    }



    selection = "within 10.0 of chain C";
    indexes = model.select_atoms(selection);
    cout << "Number of SEL atoms: " << count(indexes.begin(), indexes.end(), true) << endl;


    selection = "backbone and residue > 796 and residue < 848";
    indexes = model.select_atoms(selection);
    cout << "Number of SEL atoms: " << count(indexes.begin(), indexes.end(), true) << endl;

    Coor new_structure = structure.select_atoms("resname ALA GLY CYS and chain A B");
    cout << "New model size: " << new_structure.size() << endl;

    vector<string> seq_vec = structure.get_aa_sequences();
    cout << "The sequence is: ";
    for (const auto &seq : seq_vec) {
        cout << seq << endl;
    }
    cout << endl;

    vector<array<char, 2>> uniq_chain = structure.get_uniq_chain();
    cout << "Unique chains: ";
    for (const auto &chain : uniq_chain) {
        cout << chain[0] << " ";
    }
    cout << endl;

    model = structure.get_Models(0);
    vector<vector<string>> sec_vec = compute_SS(structure);
    cout << "The secondary structure is: ";

    for (const auto &sec : sec_vec) {
        for (const auto &s : sec) {
            cout << s << " ";
        }
        cout << endl;
    }

    Coor coor_1, coor_2;
    coor_1.read("3eam.pdb");
    coor_2.read("5bkg.pdb");

    std::vector<std::string> chain_1 = {"A"};
    std::vector<std::string> chain_2 = {"A"};
    std::vector<std::string> back_names = {"C", "N", "O", "CA"};
    std::string matrix_file = "";
    
    vector<int> test1, test2;
    auto result = get_common_atoms(coor_1, coor_2, {"A", "B", "C"}, {"A", "B", "C"});
    test1 = result.first;
    test2 = result.second;

    cout << "Number of common atoms in coor_1: " << test1.size() << endl;
    cout << "Number of common atoms in coor_2: " << test2.size() << endl;
    cout << "First common atom in coor_1: " << test1[0] << endl;
    cout << "First common atom in coor_2: " << test2[0] << endl;
    cout << "Last common atom in coor_1: " << test1[test1.size() - 1] << endl;
    cout << "Last common atom in coor_2: " << test2[test2.size() - 1] << endl;

    // string test_seq_1 = "AAGCTGAC";
    // string test_seq_2 = "AAGCTGACG";
    // Alignment align_test = sequence_align(test_seq_1, test_seq_2);
    // print_alignment(align_test);


    return 0;
}
