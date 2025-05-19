#include <cstring>
#include <iomanip>

#include "Model.h"

using namespace std;

void Model::clear() {
    x_.clear(); y_.clear(); z_.clear();
    name_.clear(); resname_.clear();
    resid_.clear(); chain_.clear();
}

size_t Model::size() const {
    return x_.size();
}

bool Model::addAtom(
           int num,
           const array<char, 5>& name_array,
           const array<char, 5>& resname_array,
           int res_id,
           const array<char, 2>& chain_array,
           float x, float y, float z, float occ, float beta,
           const array<char, 2>& alterloc = {' ', '\0'},
           const array<char, 5>& elem = {' ', ' ', ' ', ' ', '\0'},
           const array<char, 2>& insertres = {' ', '\0'},
           bool field = true) {

    field_.push_back(field);
    num_.push_back(num);
    name_.push_back(name_array);
    alterloc_.push_back(alterloc);
    resname_.push_back(resname_array);
    chain_.push_back(chain_array);
    insertres_.push_back(insertres);
    resid_.push_back(res_id);
    x_.push_back(x);
    y_.push_back(y);
    z_.push_back(z);
    occ_.push_back(occ);
    beta_.push_back(beta);
    elem_.push_back(elem);

    return true;
}