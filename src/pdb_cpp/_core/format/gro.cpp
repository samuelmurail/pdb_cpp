// GRO (GROMACS) format reader/writer for pdb_cpp
//
// GRO format reference:
//   https://manual.gromacs.org/documentation/current/reference-manual/file-formats.html#gro
//
// Fixed-width columns:
//   %5d  residue number
//   %-5s  residue name
//   %5s  atom name (right-justified)
//   %5d  atom number
//   %8.3f %8.3f %8.3f  x, y, z  (nanometres)
//   [%8.4f %8.4f %8.4f]  vx, vy, vz (optional, nm/ps)
//
// Last line of each frame: box vectors (3 or 9 values, nanometres)
//
// Coordinates are stored internally in Ångströms (1 nm = 10 Å).

#include <cctype>
#include <cstring>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>

#include "../Model.h"
#include "../Coor.h"

using namespace std;

namespace {

string gro_trim(const string& value) {
    size_t start = value.find_first_not_of(' ');
    if (start == string::npos) return "";
    size_t end = value.find_last_not_of(' ');
    return value.substr(start, end - start + 1);
}

// Derive element symbol from atom name (heuristic: first letter that isn't a digit)
string gro_derive_element(const string& atom_name) {
    if (atom_name.empty()) return "";
    // Skip leading digits
    size_t pos = 0;
    while (pos < atom_name.size() && isdigit(static_cast<unsigned char>(atom_name[pos]))) ++pos;
    if (pos >= atom_name.size()) return "";
    if (pos + 1 < atom_name.size() && islower(static_cast<unsigned char>(atom_name[pos + 1]))) {
        return atom_name.substr(pos, 2);
    }
    return string(1, atom_name[pos]);
}

} // anonymous namespace

// ---------------------------------------------------------------------------
// GRO reader
// ---------------------------------------------------------------------------

Coor GRO_parse(const string& filename) {

    Coor coor;
    coor.clear();

    ifstream file(filename);
    if (!file) {
        cerr << "Error: cannot open file " << filename << endl;
        return coor;
    }

    string line;
    int num_atoms = 0;
    int line_index = 0;      // position within current frame (0 = title, 1 = count, 2..count+1 = atoms, count+2 = box)
    int uniq_resid = -1;
    int old_resid = -99999999;

    Model model;
    model.clear();

    while (getline(file, line)) {
        if (line_index == 0) {
            // Title line — skip
            ++line_index;
            continue;
        }
        if (line_index == 1) {
            // Atom count
            num_atoms = stoi(gro_trim(line));
            ++line_index;
            continue;
        }
        if (line_index >= 2 && line_index < 2 + num_atoms) {
            // Atom line
            // Columns (0-indexed):  0-4 resid,  5-9 resname,  10-14 atom name,  15-19 atom number,
            //                       20-27 x,  28-35 y,  36-43 z   (nm)

            int resid    = stoi(line.substr(0, 5));
            string rname = gro_trim(line.substr(5, 5));
            string aname = gro_trim(line.substr(10, 5));
            int atom_num = stoi(gro_trim(line.substr(15, 5)));
            float x = stof(line.substr(20, 8)) * 10.0f;  // nm → Å
            float y = stof(line.substr(28, 8)) * 10.0f;
            float z = stof(line.substr(36, 8)) * 10.0f;

            if (resid != old_resid) {
                ++uniq_resid;
                old_resid = resid;
            }

            // Build arrays
            array<char, 5> name_array{};
            for (size_t i = 0; i < aname.size() && i < 4; ++i) name_array[i] = aname[i];
            name_array[min(aname.size(), (size_t)4)] = '\0';

            array<char, 5> resname_array{};
            for (size_t i = 0; i < rname.size() && i < 4; ++i) resname_array[i] = rname[i];
            resname_array[min(rname.size(), (size_t)4)] = '\0';

            array<char, 2> chain_array   = {'?', '\0'};  // GRO has no chain ID
            array<char, 2> alterloc      = {' ', '\0'};
            array<char, 2> insertres     = {' ', '\0'};

            // Derive element
            string elem_str = gro_derive_element(aname);
            array<char, 5> elem{};
            for (size_t i = 0; i < elem_str.size() && i < 4; ++i) elem[i] = elem_str[i];
            elem[elem_str.size()] = '\0';

            model.addAtom(
                atom_num,
                name_array,
                resname_array,
                resid,
                chain_array,
                x, y, z,
                0.0f,     // occ (not in GRO)
                0.0f,     // beta (not in GRO)
                alterloc,
                elem,
                insertres,
                false,    // field = ATOM
                uniq_resid
            );

            ++line_index;
            continue;
        }
        if (line_index == 2 + num_atoms) {
            // Box vectors line — parse into crystal_pack (orthorhombic approximation)
            // GRO box: v1x v2y v3z [v1y v1z v2x v2z v3x v3y]  all in nm
            istringstream iss(line);
            float v1x = 0, v2y = 0, v3z = 0;
            if (iss >> v1x >> v2y >> v3z) {
                // Convert nm → Å and build a CRYST1-like representation
                // For orthorhombic boxes: a=v1x, b=v2y, c=v3z, angles=90
                ostringstream cryst;
                cryst << "CRYST1"
                      << setw(9) << fixed << setprecision(3) << v1x * 10.0f
                      << setw(9) << fixed << setprecision(3) << v2y * 10.0f
                      << setw(9) << fixed << setprecision(3) << v3z * 10.0f
                      << setw(7) << fixed << setprecision(2) << 90.0f
                      << setw(7) << fixed << setprecision(2) << 90.0f
                      << setw(7) << fixed << setprecision(2) << 90.0f
                      << " P 1         1\n";
                coor.crystal_pack.set_CRYST1_pdb(cryst.str());
            }

            // Save current model
            if (model.size() > 0) {
                coor.add_Model(model);
                model.clear();
                uniq_resid = -1;
                old_resid = -99999999;
            }

            // Reset for next frame
            line_index = 0;
            continue;
        }
        ++line_index;
    }
    // In case file ends without a box line
    if (model.size() > 0) {
        coor.add_Model(model);
    }
    return coor;
}

// ---------------------------------------------------------------------------
// GRO writer
// ---------------------------------------------------------------------------

static string get_gro_string(const Coor& coor) {
    ostringstream oss;

    for (size_t model_index = 0; model_index < coor.model_size(); ++model_index) {
        const Model& model = coor.get_Models(model_index);
        oss << "Created with pdb_cpp\n";
        oss << setw(5) << model.size() << "\n";

        for (size_t i = 0; i < model.size(); ++i) {
            int resid = model.get_resid()[i];
            // GRO resid wraps at 99999
            if (resid > 99999) resid %= 100000;
            int atom_num = model.get_num()[i];
            if (atom_num > 99999) atom_num %= 100000;

            // Residue name (left-justified, 5 chars)
            string rname(model.get_resname()[i].data());
            // Atom name (right-justified, 5 chars)
            string aname(model.get_name()[i].data());

            oss << setw(5) << resid
                << left << setw(5) << rname
                << right << setw(5) << aname
                << setw(5) << atom_num
                << setw(8) << fixed << setprecision(3) << model.get_x()[i] / 10.0f  // Å → nm
                << setw(8) << fixed << setprecision(3) << model.get_y()[i] / 10.0f
                << setw(8) << fixed << setprecision(3) << model.get_z()[i] / 10.0f
                << "\n";
        }

        // Box vectors line (from crystal_pack if available, otherwise zeros)
        string cryst = coor.crystal_pack.get_pdb_crystal_pack();
        if (cryst.size() > 40) {
            // Parse back the a, b, c from the CRYST1 string
            float a = stof(cryst.substr(6, 9)) / 10.0f;   // Å → nm
            float b = stof(cryst.substr(15, 9)) / 10.0f;
            float c = stof(cryst.substr(24, 9)) / 10.0f;
            oss << setw(10) << fixed << setprecision(5) << a
                << setw(10) << fixed << setprecision(5) << b
                << setw(10) << fixed << setprecision(5) << c
                << "\n";
        } else {
            oss << "   0.00000   0.00000   0.00000\n";
        }
    }
    return oss.str();
}

bool GRO_write(const Coor& coor, const string& filename) {
    ofstream file(filename);
    if (!file) {
        cerr << "Error: cannot open file " << filename << endl;
        return false;
    }
    file << get_gro_string(coor);
    return true;
}
