#include "Coor.h"
#include "select.h"
#include <chrono>
#include <iomanip>

using namespace std;

int main() {
    Coor structure;

    auto start = chrono::high_resolution_clock::now();
    structure.read("3eam.pdb");
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

    selection = "within 10.0 of chain C";
    indexes = model.select_atoms(selection);
    cout << "Number of SEL atoms: " << count(indexes.begin(), indexes.end(), true) << endl;

    selection = "backbone and residue > 796 and residue < 848";
    indexes = model.select_atoms(selection);
    cout << "Number of SEL atoms: " << count(indexes.begin(), indexes.end(), true) << endl;

    Coor new_structure = structure.select_atoms("resname ALA GLY CYS and chain A B");
    cout << "New model size: " << new_structure.size() << endl;

    structure.get_aa_seq();

    return 0;
}
