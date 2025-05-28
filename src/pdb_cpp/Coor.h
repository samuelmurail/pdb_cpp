#ifndef COOR_H
#define COOR_H

#include <vector>
#include <string>
#include <array>

#include "Model.h"
#include "geom.h"

using namespace std;

class Coor {
public:
    // === Constructors ===
    Coor() = default;
    Coor(const string& filename) {
        read(filename);
    }
    // === Public values ===
    CrystalPack crystal_pack;
    Transformation transformation;
    Symmetry symmetry;
    size_t active_model = 0;

    // === Public interface ===

    bool read(const string& filename);
    bool write(const string & filename) const;
    //bool loadPDB(const string& filename);
    void clear();
    void add_Model(const Model& model) { models_.push_back(model); };
    size_t get_active_model() const { return active_model; };
    void set_active_model(size_t model) { active_model = model; };
    size_t size() const { return models_[active_model].size(); };
    size_t model_size() const { return models_.size(); };
    Coor select_atoms(const string &selection, size_t frame=0) const;
    Coor select_bool_index(const vector<bool> &indexes) const;
    vector<int> get_index_select(const string selection, size_t frame=0) const;
    vector<array<char, 2>> get_uniq_chain() const;
    vector<string> get_aa_sequences(bool gap_in_seq=true, size_t frame=0) const;

    Model get_Models(int i) const { return models_[i]; }
    vector<Model> get_all_Models() const { return models_; }

    // void set_crystal(float a, float b, float c, float alpha, float beta, float gamma) {
    //     crystal_pack_.set_unit_cell(alpha, beta, gamma, a, b, c);
    // }
    // CrystalPack get_crystal() const { return crystal_pack_; }


private:
    // === Storage (Structure of Arrays) ===
    vector<Model> models_;

};

#endif // COOR_H