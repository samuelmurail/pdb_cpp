#include <iostream>
#include <fstream>
#include <sstream>
#include <array>
#include <string>
#include <stdexcept>
#include <chrono>

// Constante pour la taille maximale des atomes
const size_t max_atoms = 10000; // Ajustez cette valeur en fonction de vos besoins

// Structure pour stocker les données des atomes
struct PDBData {
    std::array<std::string, max_atoms> atom_names;
    std::array<std::string, max_atoms> chains;
    std::array<std::string, max_atoms> res_names;
    std::array<double, max_atoms> x_coords;
    std::array<double, max_atoms> y_coords;
    std::array<double, max_atoms> z_coords;
    size_t count = 0; // Compteur pour le nombre d'atomes réellement lus
};

// Fonction pour parser un fichier PDB
PDBData parsePDB(const std::string& filename) {
    PDBData data;
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Unable to open file: " + filename);
    }

    std::string line;
    while (std::getline(file, line)) {
        if (line.substr(0, 4) == "ATOM" && data.count < max_atoms) {
            std::string atom_name = line.substr(12, 4);
            std::string chain = line.substr(21, 1);
            std::string res_name = line.substr(17, 3);
            double x = std::stod(line.substr(30, 8));
            double y = std::stod(line.substr(38, 8));
            double z = std::stod(line.substr(46, 8));

            data.atom_names[data.count] = atom_name;
            data.chains[data.count] = chain;
            data.res_names[data.count] = res_name;
            data.x_coords[data.count] = x;
            data.y_coords[data.count] = y;
            data.z_coords[data.count] = z;
            data.count++;
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
                  << ", Coordinates: (" << data.x_coords[i]
                  << ", " << data.y_coords[i]
                  << ", " << data.z_coords[i] << ")\n";
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

        printPDBData(data);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
