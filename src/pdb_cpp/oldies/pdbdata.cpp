#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <cstdlib> // Pour malloc et free
#include <chrono>

// Structure pour stocker les données des atomes
struct PDBData {
    std::string* atom_names;
    std::string* chains;
    std::string* res_names;
    double **xyz = new double*[3]; // Tableau 2D pour les coordonnées x, y, z
    size_t count = 0; // Compteur pour le nombre d'atomes réellement lus
};

// Fonction pour allouer de la mémoire pour les données des atomes
void allocatePDBData(PDBData& data, size_t count) {
    data.atom_names = new std::string[count];
    data.chains = new std::string[count];
    data.res_names = new std::string[count];
    for (int i = 0; i < 3; ++i) {
        data.xyz[i] = new double[count]; // Allouer de la mémoire pour les coordonnées x, y, z
    }
    data.count = count;
}

// Fonction pour libérer la mémoire allouée pour les données des atomes
void freePDBData(PDBData& data) {
    delete[] data.atom_names;
    delete[] data.chains;
    delete[] data.res_names;
    for (int i = 0; i < 3; ++i) {
        delete[] data.xyz[i]; // Free each inner array
    }
    delete[] data.xyz; // Free the outer array
    data.count = 0;
}

// Fonction pour parser un fichier PDB
PDBData parsePDB(const std::string& filename) {
    PDBData data;
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Unable to open file: " + filename);
    }

    // Compter le nombre d'atomes dans le fichier
    size_t atom_count = 0;
    std::string line;
    while (std::getline(file, line)) {
        if (line.substr(0, 4) == "ATOM" || line.substr(0, 6) == "HETATM") {
            atom_count++;
        }
    }

    // Allouer de la mémoire pour les données des atomes
    allocatePDBData(data, atom_count);

    // Revenir au début du fichier
    file.clear();
    file.seekg(0, std::ios::beg);

    // Lire les données des atomes
    size_t index = 0;
    while (std::getline(file, line)) {
        if ( (line.substr(0, 4) == "ATOM" || line.substr(0, 6) == "HETATM") && index < atom_count) {
            std::string atom_name = line.substr(12, 4);
            std::string chain = line.substr(21, 1);
            std::string res_name = line.substr(17, 3);
            double x = std::stod(line.substr(30, 8));
            double y = std::stod(line.substr(38, 8));
            double z = std::stod(line.substr(46, 8));

            data.atom_names[index] = atom_name;
            data.chains[index] = chain;
            data.res_names[index] = res_name;
            data.xyz[0][index] = x;
            data.xyz[1][index] = y;
            data.xyz[2][index] = z;
            index++;
        }
    }

    return data;
}

// Fonction pour afficher les données des atomes
void printPDBData(const PDBData& data) {
    for (size_t i = 0; i < data.count; ++i) {
        std::cout << "Atom Name: " << data.atom_names[i]
                  << ", Chain: " << data.chains[i]
                  << ", Residue Name: " << data.res_names[i]
                  << ", Coordinates : (" << data.xyz[0][i]
                  << ", " << data.xyz[1][i]
                  << ", " << data.xyz[2][i] << ")\n";
    }
}

int main() {
    try {

        std::string filename = "src/pdb_cpp/tests/input/2ol9.pdb"; // Remplacez par le nom de votre fichier PDB
        
        auto start = std::chrono::high_resolution_clock::now();
        
        PDBData data = parsePDB(filename);

        auto end = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> elapsed = end - start;
        std::cout << "parsePDB Time: " << elapsed.count() << " seconds\n";
    

        start = std::chrono::high_resolution_clock::now();
        for (size_t i = 0; i < data.count; ++i) {
            data.xyz[0][i] += 1.0; // Just an example operation
            data.xyz[1][i] += 2.0; // Just an example operation
        }
        end = std::chrono::high_resolution_clock::now();
        std::cout << "Translate Time: " << elapsed.count() << " seconds\n";


        //printPDBData(data);
        freePDBData(data); // Libérer la mémoire allouée
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
