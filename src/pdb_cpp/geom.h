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

#endif // GEOM_H