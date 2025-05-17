#include <cstring>
#include <iomanip>

#include "Model.h"
#include "format/pdb.h"

using namespace std;


bool endswith (string const &fullString, string const &ending) {
    if (fullString.length() >= ending.length()) {
        return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
    } else {
        return false;
    }
}

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
bool Model::read(const string& filename) {
    ifstream file(filename);
    if (!file) {
        cerr << "Error: cannot open file " << filename << endl;
        return false;
    }
    clear();

    if (endswith(filename, ".pdb")) {
        cout << "Reading PDB file: " << filename << endl;
        *this = PDB_parse(filename);
        return true;
    } 
    return false;

}

// bool Model::loadPDB(const string& filename) {
//     ifstream file(filename);
//     if (!file) {
//         cerr << "Error: cannot open file " << filename << endl;
//         return false;
//     }

//     clear();
//     string line;

//     while (getline(file, line)) {
//         if (line.compare(0, 6, "ATOM  ") != 0 && line.compare(0, 6, "HETATM") != 0)
//             continue;
//         try {           
//             bool field = line.compare(0, 6, "ATOM  ") ? true : false;
//             int num = stoi(line.substr(6, 5));
//             array<char, 5> name_array = {line[12], line[13], line[14], line[15], '\0'};
//             array<char, 2> alterloc = {line[16], '\0'};
//             array<char, 5> resname_array = {line[17], line[18], line[19], ' ', '\0'};
//             array<char, 2> chain_array = {line[21], '\0'};
//             int res_id            = stoi(line.substr(22, 4));
//             array<char, 2> insertres = {line[26], '\0'};
//             float x = stof(line.substr(30, 8));
//             float y = stof(line.substr(38, 8));
//             float z = stof(line.substr(46, 8));
//             float occ  = stof(line.substr(54, 6));
//             float beta = stof(line.substr(60, 6));
//             array<char, 5> elem = {line[76], line[77], ' ', ' ', '\0'};

//             field_.push_back(field);
//             num_.push_back(num);
//             name_.push_back(name_array);
//             alterloc_.push_back(alterloc);
//             resname_.push_back(resname_array);
//             chain_.push_back(chain_array);
//             insertres_.push_back(insertres);
//             resid_.push_back(res_id);
//             x_.push_back(x);
//             y_.push_back(y);
//             z_.push_back(z);
//             occ_.push_back(occ);
//             beta_.push_back(beta);
//             elem_.push_back(elem);


//         } catch (...) {
//             cerr << "Warning: skipping malformed line:\n" << line << endl;
//         }
//     }

//     return true;
// }

// bool Model::writePDB(const string& filename) const {
//     ofstream file(filename);
//     if (!file) {
//         cerr << "Error: cannot open file " << filename << endl;
//         return false;
//     }
//     // cout << "Size:" << size() << endl;

//     string field = "ATOM  "; 

//     for (size_t i = 0; i < size(); ++i) {
//         //cout << "Writing atom " << i << " "<<num_[i] << endl;
//         field = field_[i] ? "HETATM" : "ATOM  ";
//         file << field
//              << setw(6) << num_[i] << " "
//              << setw(4) << name_[i].data()
//              << setw(1) << alterloc_[i].data()
//              << setw(3) << resname_[i].data()
//              << setw(1) << chain_[i].data()
//              << setw(4) << resid_[i]
//              << setw(1) << insertres_[i].data() << "   "
//              << setw(8) << fixed << setprecision(3) << x_[i]
//              << setw(8) << fixed << setprecision(3) << y_[i]
//              << setw(8) << fixed << setprecision(3) << z_[i]
//              << setw(6) << fixed << setprecision(2) << occ_[i]
//              << setw(6) << fixed << setprecision(2) << beta_[i] << "          "
//              << elem_[i].data()
//              << "\n";
//     }

//     return true;
// }   