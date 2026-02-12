#ifndef GEOM_H
#define GEOM_H

#pragma once
#include <cstring>
#include <iomanip>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <array>
#include <vector>

#include "usalign/Kabsch.h"

class CrystalPack {

public:

    void set_CRYST1_pdb(const std::string& line) {
        // Parse CRYST1 line
        alpha = std::stof(line.substr(6, 9));
        beta = std::stof(line.substr(15, 9));
        gamma = std::stof(line.substr(24, 9));
        a = std::stof(line.substr(33, 7));
        b = std::stof(line.substr(40, 7));
        c = std::stof(line.substr(47, 7));
        sGroup = line.substr(56, 10);
        try
        {
            z = std::stoi(line.substr(67, 3));
        }
        catch(const std::exception& e)
        {
            z = 1;
        }
        
    }

    std::string get_pdb_crystal_pack() const {
        std::stringstream ss;
        ss << "CRYST1"
           << std::setw(9) << std::setprecision(3) << std::fixed << alpha
           << std::setw(9) << std::setprecision(3) << std::fixed << beta
           << std::setw(9) << std::setprecision(3) << std::fixed << gamma
           << std::setw(7) << std::setprecision(2) << std::fixed << a
           << std::setw(7) << std::setprecision(2) << std::fixed << b
           << std::setw(7) << std::setprecision(2) << std::fixed << c
           << " P" << sGroup << "  "
           << std::setw(2) << z << "\n";
        return ss.str();
    }


    void clear() {
        alpha = beta = gamma = a = b = c = nan("");
    }

private:

    float alpha=nan(""), beta=nan(""), gamma=nan(""), a=nan(""), b=nan(""), c=nan("");
    int z;
    std::string sGroup;

};


class Transformation {
public:
    void parse_pdb_transformation(const std::string& text) {
        chains.clear();
        matrix.clear();
        std::istringstream iss(text);
        std::string line;

        while (std::getline(iss, line)) {
            if (line.size() >= 42 && line.substr(34, 7) == "CHAINS:") {
                std::string chains_str = line.substr(42);
                std::istringstream chainss(chains_str);
                std::string chain;
                while (std::getline(chainss, chain, ',')) {
                    // Remove leading/trailing whitespace
                    chain.erase(chain.begin(), std::find_if(chain.begin(), chain.end(), [](unsigned char ch) { return !std::isspace(ch); }));
                    chain.erase(std::find_if(chain.rbegin(), chain.rend(), [](unsigned char ch) { return !std::isspace(ch); }).base(), chain.end());
                    if (!chain.empty())
                        chains.push_back(chain);
                }
            } else if (line.size() >= 19 && line.substr(0, 18) == "REMARK 350   BIOMT") {
                std::istringstream vals(line.substr(19));
                std::vector<float> row;
                float val;
                while (vals >> val) {
                    row.push_back(val);
                }
                matrix.push_back(row);
            }
        }
    }
    void print() const {
        std::cout << "Chains: ";
        for (const auto& chain : chains) {
            std::cout << chain << " ";
        }
        std::cout << std::endl;

        std::cout << "Matrix:" << std::endl;
        for (const auto& row : matrix) {
            for (const auto& val : row) {
                std::cout << val << " ";
            }
            std::cout << std::endl;
        }
    }

    void clear() {
        chains.clear();
        matrix.clear();
    }
    void set_chains(const std::vector<std::string> &new_chains) {
        chains = new_chains;
    }
    void add_matrix_row(const std::vector<float> &row) {
        matrix.push_back(row);
    }
    const std::vector<std::vector<float>> &get_matrix() const {
        return matrix;
    }
private:
    std::vector<std::string> chains;
    std::vector<std::vector<float>> matrix;
};

class Symmetry {
public:
    void parse_pdb_symmetry(const std::string& text) {
        matrix.clear();
        std::istringstream iss(text);
        std::string line;
        while (std::getline(iss, line)) {
            if (line.size() >= 19 && line.substr(0, 18) == "REMARK 290   SMTRY") {
                std::istringstream vals(line.substr(19));
                std::vector<float> row;
                float val;
                while (vals >> val) {
                    row.push_back(val);
                }
                matrix.push_back(row);
            }
        }
    }
    void print() const {
        std::cout << "Symmetry Matrix:" << std::endl;
        for (const auto& row : matrix) {
            for (const auto& val : row) {
                std::cout << val << " ";
            }
            std::cout << std::endl;
        }
    }
    void clear() {
        matrix.clear();
    }
private:
    std::vector<std::vector<float>> matrix;
};

inline float calculate_distance(float x1, float y1, float z1, float x2, float y2, float z2){
    return std::sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
}
inline float calculate_square_distance(float x1, float y1, float z1, float x2, float y2, float z2) {
    return (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2);
}

inline float atom_dihed_angle(
    const std::array<float, 3> &a,
    const std::array<float, 3> &b,
    const std::array<float, 3> &c,
    const std::array<float, 3> &d) {
    float ab[3] = {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
    float bc[3] = {c[0] - b[0], c[1] - b[1], c[2] - b[2]};
    float cd[3] = {d[0] - c[0], d[1] - c[1], d[2] - c[2]};

    float v1[3] = {
        ab[1] * bc[2] - ab[2] * bc[1],
        ab[2] * bc[0] - ab[0] * bc[2],
        ab[0] * bc[1] - ab[1] * bc[0]
    };
    float v2[3] = {
        cd[1] * bc[2] - cd[2] * bc[1],
        cd[2] * bc[0] - cd[0] * bc[2],
        cd[0] * bc[1] - cd[1] * bc[0]
    };
    float v1_x_v2[3] = {
        v1[1] * v2[2] - v1[2] * v2[1],
        v1[2] * v2[0] - v1[0] * v2[2],
        v1[0] * v2[1] - v1[1] * v2[0]
    };

    float bc_norm = std::sqrt(bc[0] * bc[0] + bc[1] * bc[1] + bc[2] * bc[2]);
    if (bc_norm == 0.0f) {
        return 0.0f;
    }

    float y = (v1_x_v2[0] * bc[0] + v1_x_v2[1] * bc[1] + v1_x_v2[2] * bc[2]) / bc_norm;
    float x = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
    return std::atan2(y, x) * 180.0f / static_cast<float>(M_PI);
}

inline void distance_matrix(const float* a, size_t n, const float* b, size_t m, float* out) {
    for (size_t i = 0; i < n; ++i) {
        size_t a_offset = i * 3;
        for (size_t j = 0; j < m; ++j) {
            size_t b_offset = j * 3;
            float dx = a[a_offset] - b[b_offset];
            float dy = a[a_offset + 1] - b[b_offset + 1];
            float dz = a[a_offset + 2] - b[b_offset + 2];
            out[i * m + j] = std::sqrt(dx * dx + dy * dy + dz * dz);
        }
    }
}

inline std::vector<std::vector<float>> distance_matrix(
    const std::vector<std::array<float, 3>> &a,
    const std::vector<std::array<float, 3>> &b) {
    std::vector<std::vector<float>> out(a.size(), std::vector<float>(b.size(), 0.0f));
    for (size_t i = 0; i < a.size(); ++i) {
        for (size_t j = 0; j < b.size(); ++j) {
            float dx = a[i][0] - b[j][0];
            float dy = a[i][1] - b[j][1];
            float dz = a[i][2] - b[j][2];
            out[i][j] = std::sqrt(dx * dx + dy * dy + dz * dz);
        }
    }
    return out;
}

// Wrapper function to use Kabsch algorithm with std::vector and std::array
// This maintains compatibility with the existing quaternion_rotate interface
inline std::array<std::array<float, 3>, 3> quaternion_rotate(
    const std::vector<std::array<float, 3>>& X, 
    const std::vector<std::array<float, 3>>& Y) {
    
    if (X.size() != Y.size()) {
        throw std::runtime_error("Input arrays X and Y must have the same length");
    }
    
    size_t n = X.size();
    if (n == 0) {
        // Return identity matrix for empty input
        return {{{{1.0f, 0.0f, 0.0f}}, {{0.0f, 1.0f, 0.0f}}, {{0.0f, 0.0f, 1.0f}}}};
    }
    
    // Allocate memory for Kabsch algorithm
    double **x = new double*[n];
    double **y = new double*[n];
    for (size_t i = 0; i < n; i++) {
        x[i] = new double[3];
        y[i] = new double[3];
        for (int j = 0; j < 3; j++) {
            x[i][j] = static_cast<double>(X[i][j]);
            y[i][j] = static_cast<double>(Y[i][j]);
        }
    }
    
    double rms;
    double t[3];
    double u[3][3];
    
    // Call Kabsch algorithm (mode=1 to calculate rotation matrix only)
    bool success = Kabsch(x, y, static_cast<int>(n), 1, &rms, t, u);
    
    // Clean up memory
    for (size_t i = 0; i < n; i++) {
        delete[] x[i];
        delete[] y[i];
    }
    delete[] x;
    delete[] y;
    
    if (!success) {
        // Return identity matrix if Kabsch fails
        return {{{{1.0f, 0.0f, 0.0f}}, {{0.0f, 1.0f, 0.0f}}, {{0.0f, 0.0f, 1.0f}}}};
    }
    
    // Convert double matrix to float array
    std::array<std::array<float, 3>, 3> rotation_matrix;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            rotation_matrix[i][j] = static_cast<float>(u[i][j]);
        }
    }
    
    return rotation_matrix;
}

#endif // GEOM_H