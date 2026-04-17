#ifndef COOR_H
#define COOR_H

#include <vector>
#include <string>
#include <array>
#include <unordered_map>

#include "Model.h"
#include "geom.h"

class Coor {
public:
    // === Constructors ===
    Coor() = default;
    Coor(const std::string& filename) {
        read(filename);
    }
    // === Public values ===
    CrystalPack crystal_pack;
    Transformation transformation;
    Symmetry symmetry;
    size_t active_model = 0;
    std::unordered_map<int, std::vector<int>> conect;

    // === Public interface ===

    bool read(const std::string& filename, const std::string& format = "");
    bool write(const std::string & filename) const;
    //bool loadPDB(const string& filename);
    void clear();
    void add_Model(const Model& model) { models_.push_back(model); };
    size_t get_active_model() const { return active_model; };
    void set_active_model(size_t model) { active_model = model; };
    size_t size() const {
        if (models_.empty() || active_model >= models_.size()) {
            return 0;
        }
        return models_[active_model].size();
    };
    size_t model_size() const { return models_.size(); };
    Coor select_atoms(const std::string &selection, size_t frame=0) const;
    Coor select_bool_index(const std::vector<bool> &indexes) const;
    std::vector<int> get_index_select(const std::string selection, size_t frame=0) const;
    std::vector<std::array<char, 2>> get_uniq_chain() const;
    std::vector<std::string> get_uniq_chain_str() const;
    std::vector<std::string> get_aa_sequences(bool gap_in_seq=true, size_t frame=0) const;
    std::vector<std::string> get_aa_sequences_dl(bool gap_in_seq=true, size_t frame=0) const;
    std::unordered_map<std::string, std::string> get_aa_seq(bool gap_in_seq=true, size_t frame=0) const;
    std::unordered_map<std::string, std::string> get_aa_DL_seq(bool gap_in_seq=true, size_t frame=0) const;
    std::unordered_map<std::string, std::string> get_aa_na_seq(bool gap_in_seq=true, size_t frame=0) const;
    Coor remove_incomplete_backbone_residues(const std::vector<std::string> &back_atom) const;

    Model get_Models(int i) const { return models_[i]; }
    std::vector<Model> get_all_Models() const { return models_; }

    // Per-atom getters (active model)
    const std::vector<float> &get_x(size_t frame = 0) const { return models_[frame].get_x(); }
    const std::vector<float> &get_y(size_t frame = 0) const { return models_[frame].get_y(); }
    const std::vector<float> &get_z(size_t frame = 0) const { return models_[frame].get_z(); }
    const std::vector<float> &get_beta(size_t frame = 0) const { return models_[frame].get_beta(); }
    const std::vector<float> &get_occ(size_t frame = 0) const { return models_[frame].get_occ(); }
    const std::vector<int> &get_num(size_t frame = 0) const { return models_[frame].get_num(); }
    const std::vector<int> &get_resid(size_t frame = 0) const { return models_[frame].get_resid(); }
    const std::vector<int> &get_uniqresid(size_t frame = 0) const { return models_[frame].get_uniqresid(); }
    const std::vector<std::array<char, 5>> &get_name(size_t frame = 0) const { return models_[frame].get_name(); }
    const std::vector<std::array<char, 5>> &get_resname(size_t frame = 0) const { return models_[frame].get_resname(); }
    const std::vector<std::array<char, 5>> &get_elem(size_t frame = 0) const { return models_[frame].get_elem(); }
    const std::vector<std::array<char, 2>> &get_chain(size_t frame = 0) const { return models_[frame].get_chain(); }
    const std::vector<std::array<char, 2>> &get_alterloc(size_t frame = 0) const { return models_[frame].get_alterloc(); }
    const std::vector<std::array<char, 2>> &get_insertres(size_t frame = 0) const { return models_[frame].get_insertres(); }

    // Per-atom setters (active model) — write directly, no copy
    void set_x(size_t index, float v)  { models_[active_model].set_x(index, v); }
    void set_y(size_t index, float v)  { models_[active_model].set_y(index, v); }
    void set_z(size_t index, float v)  { models_[active_model].set_z(index, v); }
    void set_beta(size_t index, float v) { models_[active_model].set_beta(index, v); }
    void set_occ(size_t index, float v)  { models_[active_model].set_occ(index, v); }
    void set_num(size_t index, int v)    { models_[active_model].set_num(index, v); }
    void set_resid(size_t index, int v)  { models_[active_model].set_resid(index, v); }
    void set_uniqresid(size_t index, int v) { models_[active_model].set_uniqresid(index, v); }
    void set_name(size_t index, const std::string &v)      { models_[active_model].set_name(index, v); }
    void set_resname(size_t index, const std::string &v)   { models_[active_model].set_resname(index, v); }
    void set_elem(size_t index, const std::string &v)      { models_[active_model].set_elem(index, v); }
    void set_chain(size_t index, const std::string &v)     { models_[active_model].set_chain(index, v); }
    void set_alterloc(size_t index, const std::string &v)  { models_[active_model].set_alterloc(index, v); }
    void set_insertres(size_t index, const std::string &v) { models_[active_model].set_insertres(index, v); }

    // void set_crystal(float a, float b, float c, float alpha, float beta, float gamma) {
    //     crystal_pack_.set_unit_cell(alpha, beta, gamma, a, b, c);
    // }
    // CrystalPack get_crystal() const { return crystal_pack_; }


private:
    // === Storage (Structure of Arrays) ===
    std::vector<Model> models_;

};

#endif // COOR_H