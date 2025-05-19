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
    int size() const { return models_.size(); }
    Model get_Models(int i) const { return models_[i]; }
    void set_crystal(float a, float b, float c, float alpha, float beta, float gamma) {
        crystal_pack_.setA(a);
        crystal_pack_.setB(b);
        crystal_pack_.setC(c);
        crystal_pack_.setAlpha(alpha);
        crystal_pack_.setBeta(beta);
        crystal_pack_.setGamma(gamma);
    }

private:
    // === Storage (Structure of Arrays) ===
    vector<Model> models_;
    CrystalPack crystal_pack_;
    
    int active_model_ = 0;

};
