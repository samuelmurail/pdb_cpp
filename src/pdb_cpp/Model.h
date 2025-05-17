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
    bool read(const string& filename);
    //bool loadPDB(const string& filename);
    bool write(const string& filename) const;
    bool addAtom(
        int num,
        const array<char, 5>& name_array,
        const array<char, 5>& resname_array,
        int res_id,
        const array<char, 2>& chain_array,
        float x, float y, float z, float occ, float beta,
        const array<char, 2>& alterloc,
        const array<char, 5>& elem,
        const array<char, 2>& insertres,
        bool field);

    void clear();
    size_t size() const;

    // Accessors
    const vector<float>& get_x() const { return x_; }
    const vector<float>& get_y() const { return y_; }
    const vector<float>& get_z() const { return z_; }
    const vector<array<char, 5>>& get_name() const { return name_; }
    const vector<array<char, 5>>& get_resname() const { return resname_; }
    const vector<int>& get_resid() const { return resid_; }
    const vector<array<char, 2>>& get_chain() const { return chain_; }
    const vector<float>& get_occ() const { return occ_; }
    const vector<float>& get_beta() const { return beta_; }
    const vector<array<char, 2>>& get_alterloc() const { return alterloc_; }
    const vector<array<char, 2>>& get_insertres() const { return insertres_; }
    const vector<array<char, 5>>& get_elem() const { return elem_; }
    const vector<int>& get_num() const { return num_; }
    const vector<bool>& get_field() const { return field_; }
    const vector<int>& get_uniqresid() const { return uniqresid_; }

private:
    // === Storage (Structure of Arrays) ===
    vector<float> x_, y_, z_, occ_, beta_;
    vector<array<char, 5>> name_, resname_, elem_;
    vector<array<char, 2>> chain_, alterloc_, insertres_;
    vector<int> num_, resid_, uniqresid_;
    vector<bool> field_;
};
