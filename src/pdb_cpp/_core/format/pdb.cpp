#include <algorithm>
#include <cctype>
#include <cstring>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>

#include "../Model.h"
#include "../Coor.h"
#include "../geom.h"
#include "encode.h"

using namespace std;

namespace {

bool is_all_spaces(const string& s) {
    for (char c : s) {
        if (c != ' ') {
            return false;
        }
    }
    return true;
}

string safe_field(const string& line, size_t start, size_t count) {
    if (start >= line.size()) {
        return string(count, ' ');
    }
    size_t len = min(count, line.size() - start);
    string out = line.substr(start, len);
    if (out.size() < count) {
        out.append(count - out.size(), ' ');
    }
    return out;
}

float parse_float_field(const string& line, size_t start, size_t count, float fallback) {
    string field = safe_field(line, start, count);
    if (is_all_spaces(field)) {
        return fallback;
    }
    return stof(field);
}

}

Coor PDB_parse(const string& filename) {

    Coor coor;
    Model model;
    string transformation_txt = "";
    string symmetry_txt = "";

    coor.clear();
    model.clear();

    ifstream file(filename);
    if (!file) {
        cerr << "Error: cannot open file " << filename << endl;
        return coor;
    }

    string line;
    int uniq_resid = -1; // Unique residue ID, will be incremented for each new residue
    int old_resid = -99999999;
    char old_insert_res = ' ';

    size_t array_i = 0;

    while (getline(file, line)) {
        if (line.compare(0, 6, "ATOM  ") == 0 || line.compare(0, 6, "HETATM") == 0) {
            bool field = line.compare(0, 6, "ATOM  ") ? true : false;
            int num = hy36decode(5, safe_field(line, 6, 5));
            array<char, 5> name_array{};
            // strip spaces
            array_i = 0;
            for (size_t i = 0; i < 4; ++i) {
                if (line[12 + i] != ' ') {
                    name_array[array_i] = line[12 + i];
                    ++array_i;
                }
            }
            name_array[array_i] = '\0';
            array<char, 2> alterloc = {line[16], '\0'};
            array<char, 5> resname_array{};
            // strip spaces
            array_i = 0;
            for (size_t i = 0; i < 3; ++i) {
                if (line[17 + i] != ' ') {
                    resname_array[array_i] = line[17 + i];
                    ++array_i;
                }
            }
            resname_array[array_i] = '\0';
            array<char, 2> chain_array = {line[21], '\0'};
            int res_id            = hy36decode(4, safe_field(line, 22, 4));
            if (res_id != old_resid || line[26] != old_insert_res) {
                ++uniq_resid;
                old_resid = res_id;
                old_insert_res = line[26];
            }
            array<char, 2> insertres = {line[26], '\0'};
            float x = parse_float_field(line, 30, 8, 0.0f);
            float y = parse_float_field(line, 38, 8, 0.0f);
            float z = parse_float_field(line, 46, 8, 0.0f);
            float occ  = parse_float_field(line, 54, 6, 1.0f);
            float beta = parse_float_field(line, 60, 6, 0.0f);
            array<char, 5> elem{};
            array_i = 0;
            for (size_t i = 0; i < 2; ++i) {
                if (line[76 + i] != ' ') {
                    elem [array_i] = line[76 + i];
                    ++array_i;
                }
            }
            elem [array_i] = '\0';


            model.addAtom(
                num,
                name_array,
                resname_array,
                res_id,
                chain_array,
                x, y, z, occ, beta,
                alterloc,
                elem,
                insertres,
                field,
                uniq_resid
            );

        } else if (line.compare(0, 3, "END") == 0) {
            if (model.size() == 0) continue;
            coor.add_Model(model);
            model.clear();
            uniq_resid = -1;
            old_resid = -99999999;
        } else if (line.compare(0, 6, "CONECT") == 0) {
            int atom_index = 0;
            bool have_atom = false;
            for (size_t pos = 6; pos + 5 <= line.size(); pos += 5) {
                string chunk = line.substr(pos, 5);
                if (is_all_spaces(chunk)) {
                    continue;
                }
                int value = hy36decode(5, chunk);
                if (!have_atom) {
                    atom_index = value;
                    have_atom = true;
                } else {
                    coor.conect[atom_index].push_back(value);
                }
            }
        } else if (line.compare(0, 6, "CRYST1") == 0) {
            // Parse CRYST1 line
            coor.crystal_pack.set_CRYST1_pdb(line);
        } else if (line.compare(0, 11, "REMARK 350 ") == 0) {
            transformation_txt += line + "\n";
        } else if (line.compare(0, 11, "REMARK 290 ") == 0) {
            symmetry_txt += line + "\n";
        }
    } 
    if (model.size() != 0) coor.add_Model(model);
    if (!transformation_txt.empty()) {
        coor.transformation.parse_pdb_transformation(transformation_txt);
        //coor.transformation = transformation;
    }
    if (!symmetry_txt.empty()) {
        coor.symmetry.parse_pdb_symmetry(symmetry_txt);
    }
    return coor;
}

string get_pdb_string(const Coor& coor) {
    ostringstream oss;
    oss << coor.crystal_pack.get_pdb_crystal_pack();

    for (size_t model_index = 0; model_index < coor.model_size(); ++model_index) {
        const Model& model = coor.get_Models(model_index);
        oss << "MODEL      " << setw(3) << model_index + 1 << "\n";
        for (size_t i = 0; i < model.size(); ++i) {
            string field = model.get_field()[i] ? "HETATM" : "ATOM  ";
            string serial = hy36encode(5, model.get_num()[i]);
            string resid = hy36encode(4, model.get_resid()[i]);
            oss << field
                << serial << " "
                << setw(4) << model.get_name()[i].data()
                << setw(1) << model.get_alterloc()[i].data()
                << setw(3) << model.get_resname()[i].data() << " "
                << setw(1) << model.get_chain()[i].data()
                << resid
                << setw(1) << model.get_insertres()[i].data() << "   "
                << setw(8) << fixed << setprecision(3) << model.get_x()[i]
                << setw(8) << fixed << setprecision(3) << model.get_y()[i]
                << setw(8) << fixed << setprecision(3) << model.get_z()[i]
                << setw(6) << fixed << setprecision(2) << model.get_occ()[i]
                << setw(6) << fixed << setprecision(2) << model.get_beta()[i] << "          "
                << model.get_elem()[i].data()
                << "\n";
        }
        oss << "ENDMDL\n";
    }
    if (!coor.conect.empty()) {
        vector<int> keys;
        keys.reserve(coor.conect.size());
        for (const auto &kv : coor.conect) {
            keys.push_back(kv.first);
        }
        sort(keys.begin(), keys.end());
        for (int atom_index : keys) {
            const vector<int> &connected_atoms = coor.conect.at(atom_index);
            for (size_t i = 0; i < connected_atoms.size(); i += 4) {
                oss << "CONECT" << hy36encode(5, atom_index);
                size_t end = min(i + 4, connected_atoms.size());
                for (size_t j = i; j < end; ++j) {
                    oss << hy36encode(5, connected_atoms[j]);
                }
                oss << "\n";
            }
        }
    }
    oss << "END\n";
    return oss.str();
}


bool PDB_write(const Coor& coor, const string& filename) {

    ofstream file(filename);
    if (!file) {
        cerr << "Error: cannot open file " << filename << endl;
        return false;
    }
    // cout << "Size:" << size() << endl;
    file << get_pdb_string(coor);

    return true;
}   