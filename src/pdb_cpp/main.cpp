#include "Coor.h"
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

    indexes = model.simple_select_atoms("resname", {"ALA", "GLY"}, "isin");
    cout << "Number of ALA GLY atoms: " << count(indexes.begin(), indexes.end(), true) << endl;
    end = chrono::high_resolution_clock::now();
    elapsed = end - start;
    cout << "Time taken to select coordinates: " << setprecision(3) << elapsed.count() << " seconds\n";

    return 0;
}
