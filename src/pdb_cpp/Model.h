#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <array>

using namespace std;

class Model {
public:
    // === Public interface ===
    bool loadPDB(const string& filename);
    void clear();
    size_t size() const;

    // Accessors
    const vector<float>& getX() const { return x_; }
    const vector<float>& getY() const { return y_; }
    const vector<float>& getZ() const { return z_; }
    const vector<array<char, 5>>& getAtomNames() const { return name_; }
    const vector<array<char, 5>>& getResNames() const { return resname_; }
    const vector<int>& getResIDs() const { return resid_; }
    const vector<array<char, 2>>& getChainIDs() const { return chain_; }

private:
    // === Storage (Structure of Arrays) ===
    vector<float> x_, y_, z_, occ_, beta_;
    vector<array<char, 5>> name_, resname_, elem_;
    vector<array<char, 2>> chain_, alterloc_, insertres_;
    vector<int> num_, resid_, uniqresid_;
    vector<bool> field_;
};
