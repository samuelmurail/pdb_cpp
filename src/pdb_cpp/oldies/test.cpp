#include <iostream>
#include <array>
#include <chrono>
#include <vector>
#include <string>

// AoS
struct Atom {
    std::string atom_name;
    std::string chain;
    std::string res_name;
    double x, y, z;
};

void benchmarkAOS() {
    const size_t max_atoms = 1000000;
    std::vector<Atom> atoms(max_atoms);

    // Remplir le vecteur
    for (size_t i = 0; i < max_atoms; ++i) {
        atoms[i] = {
            "Atom" + std::to_string(i),
            "Chain" + std::to_string(i % 2),
            "Res" + std::to_string(i % 3),
            static_cast<double>(i),
            static_cast<double>(i * 2),
            static_cast<double>(i * 3)
        };
    }

    // Mesurer le temps de parcours
    auto start = std::chrono::high_resolution_clock::now();
    double sum = 0.0;
    for (const auto& atom : atoms) {
        sum += atom.x;
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "AOS Time: " << elapsed.count() << " seconds\n";
}

// SoA
void benchmarkSOA() {
    const size_t max_atoms = 1000000;
    std::vector<std::string> atom_names(max_atoms);
    std::vector<std::string> chains(max_atoms);
    std::vector<std::string> res_names(max_atoms);
    std::vector<double> x_coords(max_atoms);
    std::vector<double> y_coords(max_atoms);
    std::vector<double> z_coords(max_atoms);

    // Remplir les vecteurs
    for (size_t i = 0; i < max_atoms; ++i) {
        atom_names[i] = "Atom" + std::to_string(i);
        chains[i] = "Chain" + std::to_string(i % 2);
        res_names[i] = "Res" + std::to_string(i % 3);
        x_coords[i] = static_cast<double>(i);
        y_coords[i] = static_cast<double>(i * 2);
        z_coords[i] = static_cast<double>(i * 3);
    }

    // Mesurer le temps de parcours
    auto start = std::chrono::high_resolution_clock::now();
    double sum = 0.0;
    for (size_t i = 0; i < max_atoms; ++i) {
        sum += x_coords[i];
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "SOA Time: " << elapsed.count() << " seconds\n";
}

int main() {
    benchmarkAOS();
    benchmarkSOA();
    return 0;
}
