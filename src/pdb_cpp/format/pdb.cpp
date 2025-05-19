#include <cstring>
#include <iomanip>

#include "../Model.h"
#include "../Coor.h"

using namespace std;

Coor PDB_parse(const string& filename) {

    Coor coor;
    Model model;

    coor.clear();
    model.clear();

    ifstream file(filename);
    if (!file) {
        cerr << "Error: cannot open file " << filename << endl;
        return coor;
    }

    string line;

    while (getline(file, line)) {
        if (line.compare(0, 6, "ATOM  ") == 0 || line.compare(0, 6, "HETATM") == 0) {
            bool field = line.compare(0, 6, "ATOM  ") ? true : false;
            int num = stoi(line.substr(6, 5));
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
                field
            );
        } else if (line.compare(0, 3, "END") == 0) {
            coor.add_Model(model);
            model.clear();
        } else if (line.compare(0, 6, "CRYST1") == 0) {
            // Parse CRYST1 line
            float a = stof(line.substr(6, 9));
            float b = stof(line.substr(15, 9));
            float c = stof(line.substr(24, 9));
            float alpha = stof(line.substr(33, 7));
            float beta = stof(line.substr(40, 7));
            float gamma = stof(line.substr(47, 7));

            coor.set_crystal(a, b, c, alpha, beta, gamma);


        }
    } 
    coor.add_Model(model);
    return coor;
}


bool PDB_write(const Coor& coor, const string& filename) {

    ofstream file(filename);
    if (!file) {
        cerr << "Error: cannot open file " << filename << endl;
        return false;
    }
    // cout << "Size:" << size() << endl;
    string field;

    for (int i = 0; i < coor.size(); ++i) {

        const Model& model = coor.get_Models(i);

        file << "MODEL     " << i + 1 << "\n";

        for (size_t i = 0; i < model.size(); ++i) {
            //cout << "Writing atom " << i << " "<<num_[i] << endl;
            field = model.get_field()[i] ? "HETATM" : "ATOM  ";
            file << field
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

        file << "ENDMDL\n";
    }

    return true;
}   