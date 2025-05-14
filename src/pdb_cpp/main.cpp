#include "Model.h"
#include <chrono>

using namespace std;

int main() {
    Model structure;

    auto start = chrono::high_resolution_clock::now();
    structure.loadPDB("src/pdb_cpp/tests/input/2ol9.pdb");

    auto end = chrono::high_resolution_clock::now();

    chrono::duration<double> elapsed = end - start;
    cout << "parsePDB Time: " << elapsed.count() << " seconds\n";

    cout << "First atom: ";
    for(int i = 0 ; i < 5 ; i ++ ){
        cout << structure.getAtomNames()[0][i] ;
    }
    for(int i = 0 ; i < 5 ; i ++ ){
        cout << structure.getResNames()[0][i] ;
    }
    // << structure.getAtomNames()[0] << " "
    // << structure.getResNames()[0] << " "
    cout << structure.getX()[0] << " "
    << structure.getY()[0] << " "
    << structure.getZ()[0] << "\n";


    return 0;
}
