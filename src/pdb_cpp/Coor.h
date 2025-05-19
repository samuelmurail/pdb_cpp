#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <array>

#include "Model.h"
#include "crystal_pack.h"

using namespace std;

class Coor {
public:
    // === Public interface ===
    bool read(const string& filename);
    bool write(const string & filename);
    //bool loadPDB(const string& filename);
    void clear();
    void add_Model(const Model& model) { models_.push_back(model); }
    Model get_Models(int i) const { return models_[i]; }

private:
    // === Storage (Structure of Arrays) ===
    vector<Model> models_;
    CrystalPack crystal_pack_;
    
    int active_model_ = 0;

};
