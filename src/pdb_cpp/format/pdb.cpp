#include <cstring>
#include <iomanip>

#include "../Model.h"
#include "../Coor.h"
#include "../geom.h"

using namespace std;

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
            int num = stoi(line.substr(6, 5));
            array<char, 5> name_array;
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
            array<char, 5> resname_array;
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
            int res_id            = stoi(line.substr(22, 4));
            if (res_id != old_resid || line[26] != old_insert_res) {
                ++uniq_resid;
                old_resid = res_id;
                old_insert_res = line[26];
            }
            array<char, 2> insertres = {line[26], '\0'};
            float x = stof(line.substr(30, 8));
            float y = stof(line.substr(38, 8));
            float z = stof(line.substr(46, 8));
            float occ  = stof(line.substr(54, 6));
            float beta = stof(line.substr(60, 6));
            array<char, 5> elem;
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
            uniq_resid = 0;
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
            std::string field = model.get_field()[i] ? "HETATM" : "ATOM  ";
            oss << field
                << setw(5) << model.get_num()[i] << " "
                << setw(4) << model.get_name()[i].data()
                << setw(1) << model.get_alterloc()[i].data()
                << setw(3) << model.get_resname()[i].data()
                << setw(1) << model.get_chain()[i].data()
                << setw(4) << model.get_resid()[i]
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