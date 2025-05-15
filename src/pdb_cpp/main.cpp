#include "Model.h"
#include <chrono>

using namespace std;

int main() {
    Model structure;

    auto start = chrono::high_resolution_clock::now();
    structure.loadPDB("3eam.pdb");

    auto end = chrono::high_resolution_clock::now();

    chrono::duration<double> elapsed = end - start;
    cout << "parsePDB Time: " << elapsed.count() << " seconds\n";

    cout << "First atom: ";
    for(int i = 0 ; i < 5 ; i ++ ){
        cout << structure.getAtomNames()[0][i] ;
    }
    cout << " ";
    for(int i = 0 ; i < 5 ; i ++ ){
        cout << structure.getResNames()[0][i] ;
    }
    cout << " ";
    // << structure.getAtomNames()[0] << " "
    // << structure.getResNames()[0] << " "
    cout << structure.getX()[0] << " "
    << structure.getY()[0] << " "
    << structure.getZ()[0] << endl;

    //structure.clear();


    start = chrono::high_resolution_clock::now();
    structure.writePDB("tmp.pdb");
    end = chrono::high_resolution_clock::now();
    elapsed = end - start;
    cout << "writePDB Time: " << elapsed.count() << " seconds\n";


    return 0;
}
