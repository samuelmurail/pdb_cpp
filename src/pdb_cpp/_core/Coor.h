#ifndef COOR_H
#define COOR_H

/**
 * @file Coor.h
 * @brief High-level structure container that owns one or more Model frames.
 */

#include <vector>
#include <string>
#include <array>
#include <unordered_map>

#include "Model.h"
#include "geom.h"

class Coor {
public:
    /** Construct an empty coordinate container. */
    Coor() = default;
    /** Load a coordinate file immediately after construction. */
    Coor(const std::string& filename) {
        read(filename);
    }

    /** Crystallographic metadata parsed from the input file. */
    CrystalPack crystal_pack;
    /** Symmetry and transformation records parsed from REMARK blocks. */
    Transformation transformation;
    /** Symmetry annotations parsed from REMARK blocks. */
    Symmetry symmetry;
    /** Index of the currently active frame. */
    size_t active_model = 0;
    /** Atom connectivity table keyed by serial number. */
    std::unordered_map<int, std::vector<int>> conect;

    /** Read a structure file into this container. */
    bool read(const std::string& filename, const std::string& format = "");
    /** Write the current structure to disk. */
    bool write(const std::string & filename) const;
    /** Remove all stored frames and metadata. */
    void clear();
    /** Append a new frame to the container. */
    void add_Model(const Model& model) { models_.push_back(model); };
    /** Return the current active frame index. */
    size_t get_active_model() const { return active_model; };
    /** Set the active frame index. */
    void set_active_model(size_t model) { active_model = model; };
    /** Return the atom count for the active frame. */
    size_t size() const {
        if (models_.empty() || active_model >= models_.size()) {
            return 0;
        }
        return models_[active_model].size();
    };
    /** Return the number of stored frames. */
    size_t model_size() const { return models_.size(); };
    /** Return a filtered copy using a selection string. */
    Coor select_atoms(const std::string &selection, size_t frame=0) const;
    /** Return a filtered copy using a boolean mask. */
    Coor select_bool_index(const std::vector<bool> &indexes) const;
    /** Return atom indices matching a selection string. */
    std::vector<int> get_index_select(const std::string selection, size_t frame=0) const;
    /** Return unique chains as fixed-size character arrays. */
    std::vector<std::array<char, 2>> get_uniq_chain() const;
    /** Return unique chains as strings. */
    std::vector<std::string> get_uniq_chain_str() const;
    /** Return amino-acid sequences per chain. */
    std::vector<std::string> get_aa_sequences(bool gap_in_seq=true, size_t frame=0) const;
    /** Return amino-acid sequences with D-residues encoded in lowercase. */
    std::vector<std::string> get_aa_sequences_dl(bool gap_in_seq=true, size_t frame=0) const;
    /** Return protein chain sequences keyed by chain identifier. */
    std::unordered_map<std::string, std::string> get_aa_seq(bool gap_in_seq=true, size_t frame=0) const;
    /** Return protein chain sequences with D-residues encoded in lowercase. */
    std::unordered_map<std::string, std::string> get_aa_DL_seq(bool gap_in_seq=true, size_t frame=0) const;
    /** Return protein and nucleic-acid sequences keyed by chain identifier. */
    std::unordered_map<std::string, std::string> get_aa_na_seq(bool gap_in_seq=true, size_t frame=0) const;
    /** Remove residues that lack the required backbone atoms. */
    Coor remove_incomplete_backbone_residues(const std::vector<std::string> &back_atom) const;

    /** Return a model by value. */
    Model get_Models(int i) const { return models_[i]; }
    /** Return all stored models by value. */
    std::vector<Model> get_all_Models() const { return models_; }

    /** Per-atom getters for the active model. */
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

    /** Per-atom setters for the active model. */
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
    /** Frame storage, each entry owning one Model. */
    std::vector<Model> models_;

};

#endif // COOR_H