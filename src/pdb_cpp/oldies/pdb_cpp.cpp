#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <chrono>

struct Atom {
    std::string recordName;
    int atomSerialNumber;
    std::string atomName;
    std::string residueName;
    char chainID;
    int residueSequenceNumber;
    float x, y, z;
};

class PDBParser {
public:
    std::vector<Atom> parse(const std::string& filePath) {
        std::vector<Atom> atoms;
        std::ifstream file(filePath);

        if (!file.is_open()) {
            std::cerr << "Error: Could not open file " << filePath << std::endl;
            return atoms;
        }

        std::string line;
        while (std::getline(file, line)) {
            if (line.substr(0, 4) == "ATOM" || line.substr(0, 6) == "HETATM") {
                Atom atom = parseAtomLine(line);
                atoms.push_back(atom);
            }
        }

        file.close();
        return atoms;
    }

private:
    Atom parseAtomLine(const std::string& line) {
        Atom atom;
        atom.recordName = line.substr(0, 6);
        atom.atomSerialNumber = std::stoi(line.substr(6, 5));
        atom.atomName = line.substr(12, 4);
        atom.residueName = line.substr(17, 3);
        atom.chainID = line[21];
        atom.residueSequenceNumber = std::stoi(line.substr(22, 4));
        atom.x = std::stof(line.substr(30, 8));
        atom.y = std::stof(line.substr(38, 8));
        atom.z = std::stof(line.substr(46, 8));
        return atom;
    }
};

int main() {
    PDBParser parser;
    std::string filePath = "src/pdb_cpp/tests/input/2ol9.pdb"; // Replace with your PDB file path

    auto start = std::chrono::high_resolution_clock::now();

    std::vector<Atom> atoms = parser.parse(filePath);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "parsePDB Time: " << elapsed.count() << " seconds\n";

    start = std::chrono::high_resolution_clock::now();
    for (auto& atom : atoms) {
        atom.x += 1.0; // Just an example operation
        atom.y += 2.0; // Just an example operation
    }
    end = std::chrono::high_resolution_clock::now();
    std::cout << "Translate Time: " << elapsed.count() << " seconds\n";

    
    for (const auto& atom : atoms) {
        std::cout << "Atom: " << atom.atomName
                  << ", Residue: " << atom.residueName
                  << ", Chain: " << atom.chainID
                  << ", Coordinates: (" << atom.x << ", " << atom.y << ", " << atom.z << ")"
                  << std::endl;
    }

    return 0;
}