#ifndef GEOM_H
#define GEOM_H

#pragma once
#include <cstring>
#include <iomanip>
#include <string>
#include <cmath>
#include <algorithm>
#include <iostream>

using namespace std;

class CrystalPack {

public:

    void set_CRYST1_pdb(const string& line) {
        // Parse CRYST1 line
        alpha = stof(line.substr(6, 9));
        beta = stof(line.substr(15, 9));
        gamma = stof(line.substr(24, 9));
        a = stof(line.substr(33, 7));
        b = stof(line.substr(40, 7));
        c = stof(line.substr(47, 7));
        sGroup = line.substr(56, 10);
        try
        {
            z = stoi(line.substr(67, 3));
        }
        catch(const exception& e)
        {
            z = 1;
        }
        
    }

    string get_pdb_crystal_pack() const {
        stringstream ss;
        ss << "CRYST1"
           << setw(9) << setprecision(3) << fixed << alpha
           << setw(9) << setprecision(3) << fixed << beta
           << setw(9) << setprecision(3) << fixed << gamma
           << setw(7) << setprecision(2) << fixed << a
           << setw(7) << setprecision(2) << fixed << b
           << setw(7) << setprecision(2) << fixed << c
           << " P" << sGroup << "  "
           << setw(2) << z << "\n";
        return ss.str();
    }


    void clear() {
        alpha = beta = gamma = a = b = c = nan("");
    }

private:

    float alpha=nan(""), beta=nan(""), gamma=nan(""), a=nan(""), b=nan(""), c=nan("");
    int z;
    string sGroup;

};


class Transformation {
public:
    void parse_pdb_transformation(const string& text) {
        chains.clear();
        matrix.clear();
        istringstream iss(text);
        string line;

        // cout << "Parsing transformation matrix..." << endl;
        // cout << "Text: " << text << endl;

        while (getline(iss, line)) {
            // cout<< "Line: " << line << endl;
            if (line.size() >= 42 && line.substr(34, 7) == "CHAINS:") {
                string chains_str = line.substr(42);
                istringstream chainss(chains_str);
                string chain;
                while (getline(chainss, chain, ',')) {
                    // Remove leading/trailing whitespace
                    chain.erase(chain.begin(), find_if(chain.begin(), chain.end(), [](unsigned char ch) { return !isspace(ch); }));
                    chain.erase(find_if(chain.rbegin(), chain.rend(), [](unsigned char ch) { return !isspace(ch); }).base(), chain.end());
                    if (!chain.empty())
                        chains.push_back(chain);
                }
            } else if (line.size() >= 19 && line.substr(0, 18) == "REMARK 350   BIOMT") {
                istringstream vals(line.substr(19));
                vector<float> row;
                float val;
                while (vals >> val) {
                    row.push_back(val);
                }
                matrix.push_back(row);
            }
        }
    }
    void print() const {
        cout << "Chains: ";
        for (const auto& chain : chains) {
            cout << chain << " ";
        }
        cout << endl;

        cout << "Matrix:" << endl;
        for (const auto& row : matrix) {
            for (const auto& val : row) {
                cout << val << " ";
            }
            cout << endl;
        }
    }

    void clear() {
        chains.clear();
        matrix.clear();
    }
private:
    vector<string> chains;
    vector<vector<float>> matrix;
};

class Symmetry {
public:
    void parse_pdb_symmetry(const string& text) {
        matrix.clear();
        istringstream iss(text);
        string line;
        while (getline(iss, line)) {
            if (line.size() >= 19 && line.substr(0, 18) == "REMARK 290   SMTRY") {
                istringstream vals(line.substr(19));
                vector<float> row;
                float val;
                while (vals >> val) {
                    row.push_back(val);
                }
                matrix.push_back(row);
            }
        }
    }
    void print() const {
        cout << "Symmetry Matrix:" << endl;
        for (const auto& row : matrix) {
            for (const auto& val : row) {
                cout << val << " ";
            }
            cout << endl;
        }
    }
    void clear() {
        matrix.clear();
    }
private:
    vector<vector<float>> matrix;
};

inline float calculate_distance(float x1, float y1, float z1, float x2, float y2, float z2){
    return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
}
inline float calculate_square_distance(float x1, float y1, float z1, float x2, float y2, float z2) {
    return (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2);
}

#endif // GEOM_H