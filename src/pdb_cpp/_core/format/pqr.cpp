// PQR format reader/writer for pdb_cpp
//
// PQR is a PDB-like format used by electrostatics tools (APBS, PDB2PQR).
// Key differences from PDB:
//   - columns 54-62 hold the partial charge (instead of occupancy)
//   - columns 62-70 hold the atomic radius  (instead of B-factor)
//   - residue name occupies columns 16-20 (4 chars, no altloc)
//   - element symbol is absent
//
// References:
//   https://pdb2pqr.readthedocs.io/en/latest/formats/pqr.html

#include <algorithm>
#include <cctype>
#include <cstring>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>

#include "../Model.h"
#include "../Coor.h"
#include "encode.h"

using namespace std;

namespace {

string pqr_safe_field(const string& line, size_t start, size_t count) {
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

bool pqr_is_all_spaces(const string& s) {
    for (char c : s) {
        if (c != ' ') return false;
    }
    return true;
}

string pqr_trim(const string& value) {
    size_t start = value.find_first_not_of(' ');
    if (start == string::npos) return "";
    size_t end = value.find_last_not_of(' ');
    return value.substr(start, end - start + 1);
}

char pqr_safe_char(const string& line, size_t index, char fallback = ' ') {
    if (index >= line.size()) return fallback;
    return line[index];
}

float pqr_parse_float(const string& line, size_t start, size_t count, float fallback) {
    string field = pqr_safe_field(line, start, count);
    if (pqr_is_all_spaces(field)) return fallback;
    return stof(field);
}

// Derive a 1- or 2-character element symbol from the atom name.
string pqr_derive_element(const string& atom_name) {
    if (atom_name.empty()) return "";
    if (isdigit(static_cast<unsigned char>(atom_name[0]))) {
        return (atom_name.size() > 1) ? string(1, atom_name[1]) : "";
    }
    if (atom_name.size() >= 2 && islower(static_cast<unsigned char>(atom_name[1]))) {
        return atom_name.substr(0, 2);
    }
    return atom_name.substr(0, 1);
}

} // anonymous namespace

// ---------------------------------------------------------------------------
// PQR reader
// ---------------------------------------------------------------------------

Coor PQR_parse(const string& filename) {

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
    int uniq_resid = -1;
    int old_resid  = -99999999;

    while (getline(file, line)) {
        if (line.compare(0, 6, "ATOM  ") == 0 || line.compare(0, 6, "HETATM") == 0) {
            bool field = line.compare(0, 6, "ATOM  ") ? true : false;

            int num = hy36decode(5, pqr_safe_field(line, 6, 5));

            // Atom name: columns 12-16
            const string name_field = pqr_safe_field(line, 12, 4);
            array<char, 5> name_array{};
            size_t ai = 0;
            for (size_t i = 0; i < 4; ++i) {
                if (name_field[i] != ' ') {
                    name_array[ai++] = name_field[i];
                }
            }
            name_array[ai] = '\0';

            // No altloc in PQR
            array<char, 2> alterloc = {' ', '\0'};

            // Residue name occupies columns 16-20 (4 chars, trimmed to 3-4)
            const string resname_field = pqr_safe_field(line, 16, 4);
            string resname_trimmed = pqr_trim(resname_field);
            // Take up to 4 characters
            array<char, 5> resname_array{};
            ai = 0;
            for (size_t i = 0; i < resname_trimmed.size() && i < 4; ++i) {
                resname_array[ai++] = resname_trimmed[i];
            }
            resname_array[ai] = '\0';

            array<char, 2> chain_array = {pqr_safe_char(line, 21), '\0'};
            int res_id = hy36decode(4, pqr_safe_field(line, 22, 4));

            if (res_id != old_resid) {
                ++uniq_resid;
                old_resid = res_id;
            }

            array<char, 2> insertres = {' ', '\0'};
            float x = pqr_parse_float(line, 30, 8, 0.0f);
            float y = pqr_parse_float(line, 38, 8, 0.0f);
            float z = pqr_parse_float(line, 46, 8, 0.0f);

            // PQR: charge in 54-62, radius in 62-70
            float charge = pqr_parse_float(line, 54, 8, 0.0f);
            float radius = pqr_parse_float(line, 62, 8, 0.0f);

            // Derive element from atom name (PQR has no element column)
            string atom_name_str = pqr_trim(name_field);
            string elem_str = pqr_derive_element(atom_name_str);
            array<char, 5> elem{};
            for (size_t i = 0; i < elem_str.size() && i < 4; ++i) {
                elem[i] = elem_str[i];
            }
            elem[elem_str.size()] = '\0';

            // Store charge in occ and radius in beta (matching pdb_numpy convention)
            model.addAtom(
                num,
                name_array,
                resname_array,
                res_id,
                chain_array,
                x, y, z,
                charge,   // occ slot ← charge
                radius,   // beta slot ← radius
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
            old_resid  = -99999999;
        } else if (line.compare(0, 6, "CRYST1") == 0) {
            coor.crystal_pack.set_CRYST1_pdb(line);
        }
    }
    if (model.size() != 0) coor.add_Model(model);
    return coor;
}

// ---------------------------------------------------------------------------
// PQR writer
// ---------------------------------------------------------------------------

static string get_pqr_string(const Coor& coor) {
    ostringstream oss;
    oss << coor.crystal_pack.get_pdb_crystal_pack();

    for (size_t model_index = 0; model_index < coor.model_size(); ++model_index) {
        const Model& model = coor.get_Models(model_index);
        oss << "MODEL      " << setw(3) << model_index + 1 << "\n";

        for (size_t i = 0; i < model.size(); ++i) {
            string field = model.get_field()[i] ? "HETATM" : "ATOM  ";
            string serial = hy36encode(5, model.get_num()[i]);
            string resid  = hy36encode(4, model.get_resid()[i]);

            // Format atom name: left-pad short names with a space
            string raw_name(model.get_name()[i].data());
            string atom_name;
            if (raw_name.size() <= 3 &&
                !raw_name.empty() &&
                (raw_name[0] == 'C' || raw_name[0] == 'H' ||
                 raw_name[0] == 'O' || raw_name[0] == 'N' ||
                 raw_name[0] == 'S' || raw_name[0] == 'P')) {
                atom_name = " " + raw_name;
            } else {
                atom_name = raw_name;
            }
            while (atom_name.size() < 4) atom_name += ' ';
            if (atom_name.size() > 4) atom_name = atom_name.substr(0, 4);

            oss << field
                << serial << " "
                << atom_name << " "
                << setw(3) << model.get_resname()[i].data() << " "
                << setw(1) << model.get_chain()[i].data()
                << resid << "    "
                << setw(8) << fixed << setprecision(3) << model.get_x()[i]
                << setw(8) << fixed << setprecision(3) << model.get_y()[i]
                << setw(8) << fixed << setprecision(3) << model.get_z()[i]
                << setw(8) << fixed << setprecision(4) << model.get_occ()[i]    // charge
                << setw(8) << fixed << setprecision(4) << model.get_beta()[i]   // radius
                << "\n";
        }
        oss << "ENDMDL\n";
    }
    oss << "END\n";
    return oss.str();
}

bool PQR_write(const Coor& coor, const string& filename) {
    ofstream file(filename);
    if (!file) {
        cerr << "Error: cannot open file " << filename << endl;
        return false;
    }
    file << get_pqr_string(coor);
    return true;
}
