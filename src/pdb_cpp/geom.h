#ifndef GEOM_H
#define GEOM_H

#pragma once
#include <cstring>
#include <iomanip>
#include <string>
#include <cmath>
#include <algorithm>
#include <iostream>

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

// makeQ function for quaternion rotation matrix
// Source: https://github.com/charnley/rmsd/blob/master/rmsd/calculate_rmsd.py
// matrix involved in quaternion rotation
inline std::array<std::array<float, 4>, 4> makeQ(float r1, float r2, float r3, float r4 = 0.0f) {
    return {{
        {{ r4, -r3,  r2,  r1}},
        {{ r3,  r4, -r1,  r2}},
        {{-r2,  r1,  r4,  r3}},
        {{-r1, -r2, -r3,  r4}}
    }};
}

// makeW function for quaternion rotation matrix
// Source: https://github.com/charnley/rmsd/blob/master/rmsd/calculate_rmsd.py
// matrix involved in quaternion rotation
inline std::array<std::array<float, 4>, 4> makeW(float r1, float r2, float r3, float r4 = 0.0f) {
    return {{
        {{ r4,  r3, -r2,  r1}},
        {{-r3,  r4,  r1,  r2}},
        {{ r2, -r1,  r4,  r3}},
        {{-r1, -r2, -r3,  r4}}
    }};
}

// Matrix transpose helper for 4x4 matrices
inline std::array<std::array<float, 4>, 4> matrix_transpose_4x4(
    const std::array<std::array<float, 4>, 4>& matrix) {
    std::array<std::array<float, 4>, 4> result{};
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            result[i][j] = matrix[j][i];
        }
    }
    return result;
}

// Matrix multiplication helper for 4x4 matrices
inline std::array<std::array<float, 4>, 4> matrix_multiply_4x4(
    const std::array<std::array<float, 4>, 4>& A, 
    const std::array<std::array<float, 4>, 4>& B) {
    std::array<std::array<float, 4>, 4> result{};
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            result[i][j] = 0.0f;
            for (int k = 0; k < 4; k++) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return result;
}

// Matrix addition helper for 4x4 matrices
inline std::array<std::array<float, 4>, 4> matrix_add_4x4(
    const std::array<std::array<float, 4>, 4>& A, 
    const std::array<std::array<float, 4>, 4>& B) {
    std::array<std::array<float, 4>, 4> result{};
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            result[i][j] = A[i][j] + B[i][j];
        }
    }
    return result;
}

// Vector dot product for 4D vectors
inline float vector_dot_4d(const std::array<float, 4>& a, const std::array<float, 4>& b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3];
}

// Vector normalization for 4D vectors
inline std::array<float, 4> vector_normalize_4d(const std::array<float, 4>& v) {
    float norm = std::sqrt(vector_dot_4d(v, v));
    if (norm < 1e-10f) {
        return {{1.0f, 0.0f, 0.0f, 0.0f}};  // Return unit vector if input is zero
    }
    return {{v[0]/norm, v[1]/norm, v[2]/norm, v[3]/norm}};
}

// Matrix-vector multiplication for 4x4 matrix and 4D vector
inline std::array<float, 4> matrix_vector_multiply_4x4(
    const std::array<std::array<float, 4>, 4>& A, 
    const std::array<float, 4>& v) {
    std::array<float, 4> result{};
    for (int i = 0; i < 4; i++) {
        result[i] = 0.0f;
        for (int j = 0; j < 4; j++) {
            result[i] += A[i][j] * v[j];
        }
    }
    return result;
}

// Power method to find largest eigenvalue and eigenvector of 4x4 symmetric matrix
// Returns pair of (eigenvalue, eigenvector)
inline std::pair<float, std::array<float, 4>> power_method_4x4(
    const std::array<std::array<float, 4>, 4>& A, 
    int max_iterations = 100, float tolerance = 1e-6f) {
    
    // Initial guess for eigenvector
    std::array<float, 4> v = {{1.0f, 1.0f, 1.0f, 1.0f}};
    v = vector_normalize_4d(v);
    
    float eigenvalue = 0.0f;
    
    for (int iter = 0; iter < max_iterations; iter++) {
        // v_new = A * v
        auto v_new = matrix_vector_multiply_4x4(A, v);
        
        // eigenvalue = v^T * A * v = v^T * v_new
        float new_eigenvalue = vector_dot_4d(v, v_new);
        
        // Normalize v_new
        v_new = vector_normalize_4d(v_new);
        
        // Check convergence
        if (std::abs(new_eigenvalue - eigenvalue) < tolerance) {
            eigenvalue = new_eigenvalue;
            v = v_new;
            break;
        }
        
        eigenvalue = new_eigenvalue;
        v = v_new;
    }
    
    return std::make_pair(eigenvalue, v);
}

// quaternion_transform function
// Source: https://github.com/charnley/rmsd/blob/master/rmsd/calculate_rmsd.py
// Get optimal rotation
// note: translation will be zero when the centroids of each molecule are the same.
inline std::array<std::array<float, 3>, 3> quaternion_transform(float r1, float r2, float r3, float r4 = 0.0f) {
    auto Wt_r = matrix_transpose_4x4(makeW(r1, r2, r3, r4));
    auto Q_r = makeQ(r1, r2, r3, r4);
    auto rot_4x4 = matrix_multiply_4x4(Wt_r, Q_r);
    
    // Extract the 3x3 rotation matrix from the 4x4 result (equivalent to [:3, :3])
    std::array<std::array<float, 3>, 3> rot{};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            rot[i][j] = rot_4x4[i][j];
        }
    }
    return rot;
}

// quaternion_rotate function
// Source: https://github.com/charnley/rmsd/blob/master/rmsd/calculate_rmsd.py
// Calculate the optimal rotation matrix between two sets of points using quaternion method
inline std::array<std::array<float, 3>, 3> quaternion_rotate(
    const std::vector<std::array<float, 3>>& X, 
    const std::vector<std::array<float, 3>>& Y) {
    
    if (X.size() != Y.size()) {
        throw std::runtime_error("Input arrays X and Y must have the same length");
    }
    
    size_t N = X.size();
    if (N == 0) {
        // Return identity matrix for empty input
        return {{{{1.0f, 0.0f, 0.0f}}, {{0.0f, 1.0f, 0.0f}}, {{0.0f, 0.0f, 1.0f}}}};
    }
    
    // Create W matrices for each point in Y
    std::vector<std::array<std::array<float, 4>, 4>> W_matrices;
    W_matrices.reserve(N);
    for (size_t k = 0; k < N; k++) {
        W_matrices.push_back(makeW(Y[k][0], Y[k][1], Y[k][2]));
    }
    
    // Create Q matrices for each point in X
    std::vector<std::array<std::array<float, 4>, 4>> Q_matrices;
    Q_matrices.reserve(N);
    for (size_t k = 0; k < N; k++) {
        Q_matrices.push_back(makeQ(X[k][0], X[k][1], X[k][2]));
    }
    
    // Compute Qt_dot_W for each point (Q[k].T * W[k])
    std::vector<std::array<std::array<float, 4>, 4>> Qt_dot_W_matrices;
    Qt_dot_W_matrices.reserve(N);
    for (size_t k = 0; k < N; k++) {
        auto Qt = matrix_transpose_4x4(Q_matrices[k]);
        Qt_dot_W_matrices.push_back(matrix_multiply_4x4(Qt, W_matrices[k]));
    }
    
    // Sum all Qt_dot_W matrices to create matrix A
    std::array<std::array<float, 4>, 4> A{};
    // Initialize A to zero
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            A[i][j] = 0.0f;
        }
    }
    
    // Sum all matrices
    for (size_t k = 0; k < N; k++) {
        A = matrix_add_4x4(A, Qt_dot_W_matrices[k]);
    }
    
    // Find largest eigenvalue and corresponding eigenvector
    auto eigen_result = power_method_4x4(A);
    std::array<float, 4> r = eigen_result.second;
    
    // Convert quaternion to rotation matrix using quaternion_transform
    return quaternion_transform(r[0], r[1], r[2], r[3]);
}

#endif // GEOM_H