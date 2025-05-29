#ifndef MODEL_H
#define MODEL_H

#pragma once
#include <vector>
#include <string>
#include <array>

struct Token;

class Model
{
public:
    // === Public interface ===
    // bool loadPDB(const std::string& filename);
    bool addAtom(
        int num,
        const std::array<char, 5> &name_array,
        const std::array<char, 5> &resname_array,
        int res_id,
        const std::array<char, 2> &chain_array,
        float x, float y, float z, float occ, float beta,
        const std::array<char, 2> &alterloc,
        const std::array<char, 5> &elem,
        const std::array<char, 2> &insertres,
        bool field,
        int uniqresid);

    void clear();
    size_t size() const;
    std::vector<bool> simple_select_atoms(const std::string &column, const std::vector<std::string> &values, const std::string &op) const;
    std::vector<bool> select_tokens(const Token &tokens) const;
    std::vector<bool> select_atoms(const std::string selection) const;
    Model select_index(const std::vector<bool> &indexes) const;
    std::vector<int> get_index_select(const std::string selection) const;

    // Accessors
    const std::vector<float> &get_x() const { return x_; }
    const std::vector<float> &get_y() const { return y_; }
    const std::vector<float> &get_z() const { return z_; }
    const std::vector<std::array<char, 5>> &get_name() const { return name_; }
    const std::vector<std::array<char, 5>> &get_resname() const { return resname_; }
    const std::vector<int> &get_resid() const { return resid_; }
    const std::vector<std::array<char, 2>> &get_chain() const { return chain_; }
    const std::vector<float> &get_occ() const { return occ_; }
    const std::vector<float> &get_beta() const { return beta_; }
    const std::vector<std::array<char, 2>> &get_alterloc() const { return alterloc_; }
    const std::vector<std::array<char, 2>> &get_insertres() const { return insertres_; }
    const std::vector<std::array<char, 5>> &get_elem() const { return elem_; }
    const std::vector<int> &get_num() const { return num_; }
    const std::vector<bool> &get_field() const { return field_; }
    const std::vector<int> &get_uniqresid() const { return uniqresid_; }

    std::vector<std::array<char, 2>> get_uniq_chain() const;

    float distance(size_t i, size_t j) const;
    
    // Calculate centroid of all atoms or specific indices
    std::array<float, 3> get_centroid() const;
    std::array<float, 3> get_centroid(const std::vector<int>& indices) const;

    // Coordinate setters
    void set_x(size_t index, float value) { if (index < x_.size()) x_[index] = value; }
    void set_y(size_t index, float value) { if (index < y_.size()) y_[index] = value; }
    void set_z(size_t index, float value) { if (index < z_.size()) z_[index] = value; }

private:
    // === Storage (Structure of Arrays) ===
    std::vector<float> x_, y_, z_, occ_, beta_;
    std::vector<std::array<char, 5>> name_, resname_, elem_;
    std::vector<std::array<char, 2>> chain_, alterloc_, insertres_;
    std::vector<int> num_, resid_, uniqresid_;
    std::vector<bool> field_;
};

#endif // MODEL_H