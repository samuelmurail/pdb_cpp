#include <cstring>

#include "Model.h"

using namespace std;

void Model::clear() {
    x_.clear(); y_.clear(); z_.clear();
    name_.clear(); resname_.clear();
    resid_.clear(); chain_id_.clear();
}

size_t Model::size() const {
    return x_.size();
}

void Model::trimSpaces(string& s) {
    s.erase(remove_if(s.begin(), s.end(), ::isspace), s.end());
}

bool Model::loadPDB(const string& filename) {
    ifstream file(filename);
    if (!file) {
        cerr << "Error: cannot open file " << filename << "\n";
        return false;
    }

    clear();
    string line;

    while (getline(file, line)) {
        if (line.compare(0, 6, "ATOM  ") != 0 && line.compare(0, 6, "HETATM") != 0)
            continue;

        try {
            char chain_id         = line[21];
            int res_id            = stoi(line.substr(22, 4));

            float x = stof(line.substr(30, 8));
            float y = stof(line.substr(38, 8));
            float z = stof(line.substr(46, 8));


            x_.push_back(x);
            y_.push_back(y);
            z_.push_back(z);

            array<char, 5> resname_array = {' ', ' ', ' ', ' ', '\0'};
            memcpy(resname_array.data(), line.substr(17, 3).c_str(), 4);
            
            array<char, 5> name_array = {' ', ' ', ' ', ' ', '\0'};
            memcpy(name_array.data(), line.substr(12, 4).c_str(), 4);
            

            name_.push_back(name_array);
            resname_.push_back(resname_array);
            resid_.push_back(res_id);
            chain_id_.push_back(chain_id);

        } catch (...) {
            cerr << "Warning: skipping malformed line:\n" << line << "\n";
        }
    }

    return true;
}