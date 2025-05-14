#pragma once
#include <vector>
#include <string>

#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <chrono>

struct AtomicCoordinates {
    std::vector<float> x, y, z;            // Cartesian coordinates
    std::vector<std::string> atom_name;    // Atom names (e.g. "CA", "N")
    std::vector<std::string> res_name;     // Residue names (e.g. "GLY", "ALA")
    std::vector<int> res_id;               // Residue sequence numbers
    std::vector<char> chain_id;            // Chain ID (e.g. 'A', 'B')

    size_t size() const { return x.size(); }

    void clear() {
        x.clear(); y.clear(); z.clear();
        atom_name.clear(); res_name.clear();
        res_id.clear(); chain_id.clear();
    }
};



bool parsePDB(const std::string& filename, AtomicCoordinates& coords) {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Cannot open file: " << filename << "\n";
        return false;
    }

    coords.clear();
    std::string line;

    while (std::getline(file, line)) {
        if (line.compare(0, 6, "ATOM  ") != 0 && line.compare(0, 6, "HETATM") != 0)
            continue;

        try {
            std::string atom_name = line.substr(12, 4);
            std::string res_name  = line.substr(17, 3);
            char chain_id         = line[21];
            int res_id            = std::stoi(line.substr(22, 4));

            float x = std::stof(line.substr(30, 8));
            float y = std::stof(line.substr(38, 8));
            float z = std::stof(line.substr(46, 8));

            // Trim spaces from atom/residue names
            atom_name.erase(remove_if(atom_name.begin(), atom_name.end(), isspace), atom_name.end());
            res_name.erase(remove_if(res_name.begin(), res_name.end(), isspace), res_name.end());

            coords.x.push_back(x);
            coords.y.push_back(y);
            coords.z.push_back(z);
            coords.atom_name.push_back(atom_name);
            coords.res_name.push_back(res_name);
            coords.res_id.push_back(res_id);
            coords.chain_id.push_back(chain_id);

        } catch (...) {
            std::cerr << "Failed to parse line: " << line << "\n";
        }
    }

    return true;
}

int main() {
    AtomicCoordinates atoms;

    auto start = std::chrono::high_resolution_clock::now();
    parsePDB("src/pdb_cpp/tests/input/2ol9.pdb", atoms);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "parsePDB Time: " << elapsed.count() << " seconds\n";

}