#include <cstring>

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

bool Model::loadPDB(const string& filename) {
    ifstream file(filename);
    if (!file) {
        cerr << "Error: cannot open file " << filename << endl;
        return false;
    }

    clear();
    string line;

    while (getline(file, line)) {
        if (line.compare(0, 6, "ATOM  ") != 0 && line.compare(0, 6, "HETATM") != 0)
            continue;
        try {           
            bool field = line.compare(0, 6, "ATOM  ") ? true : false;
            int num = stoi(line.substr(7, 5));
            array<char, 5> name_array = {line[12], line[13], line[14], line[15], '\0'};
            array<char, 2> alterloc = {line[16], '\0'};
            array<char, 5> resname_array = {line[17], line[18], line[19], ' ', '\0'};
            array<char, 2> chain_array = {line[21], '\0'};
            int res_id            = stoi(line.substr(22, 4));
            array<char, 2> insertres = {line[26], '\0'};
            float x = stof(line.substr(30, 8));
            float y = stof(line.substr(38, 8));
            float z = stof(line.substr(46, 8));
            float occ  = stof(line.substr(54, 6));
            float beta = stof(line.substr(60, 6));
            array<char, 5> elem = {line[76], line[77], ' ', ' ', '\0'};

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


        } catch (...) {
            cerr << "Warning: skipping malformed line:\n" << line << endl;
        }
    }

    return true;
}