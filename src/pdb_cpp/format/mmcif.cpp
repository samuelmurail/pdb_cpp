#include <algorithm>
#include <cctype>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "mmcif.h"
#include "../Coor.h"
#include "../Model.h"
#include "../geom.h"

using namespace std;

namespace {

bool starts_with(const string &line, const string &prefix) {
    return line.size() >= prefix.size() && line.compare(0, prefix.size(), prefix) == 0;
}

string trim_quotes(const string &value) {
    if (value.size() >= 2) {
        if ((value.front() == '"' && value.back() == '"') ||
            (value.front() == '\'' && value.back() == '\'')) {
            return value.substr(1, value.size() - 2);
        }
    }
    return value;
}

vector<string> split_tokens(const string &line) {
    vector<string> tokens;
    string token;
    bool in_quote = false;
    char quote_char = '\0';

    for (size_t i = 0; i < line.size(); ++i) {
        char ch = line[i];
            // Stop at inline comments when not inside quotes.
            if (!in_quote && ch == '#') {
            break;
        }
        if (in_quote) {
            if (ch == quote_char) {
                in_quote = false;
                tokens.push_back(token);
                token.clear();
            } else {
                token += ch;
            }
            continue;
        }
        // Keep quoted tokens together, excluding the quote characters.
        if (ch == '\'' || ch == '"') {
            in_quote = true;
            quote_char = ch;
            continue;
        }
        if (isspace(static_cast<unsigned char>(ch))) {
            if (!token.empty()) {
                tokens.push_back(token);
                token.clear();
            }
            continue;
        }
        token += ch;
    }
    if (!token.empty()) {
        tokens.push_back(token);
    }
    return tokens;
}

struct LoopTable {
    vector<string> col_names;
    vector<vector<string>> rows;
};

LoopTable parse_loop_for_prefix(const vector<string> &lines, const string &prefix) {
    LoopTable table;
    bool in_loop = false;
    vector<string> pending;

    // Track multi-line tokens wrapped by ';' delimiters.
    bool in_multiline = false;
    string multiline_value;

    for (size_t i = 0; i < lines.size(); ++i) {
        string line = lines[i];
        if (!line.empty() && line.back() == '\n') {
            line.pop_back();
        }

        if (in_multiline) {
            if (starts_with(line, ";")) {
                in_multiline = false;
                string token = multiline_value;
                multiline_value.clear();
                pending.push_back(token);
                if (table.col_names.size() > 0 && pending.size() == table.col_names.size()) {
                    table.rows.push_back(pending);
                    pending.clear();
                }
                continue;
            }
            multiline_value += line + "\n";
            continue;
        }

        if (starts_with(line, "loop_")) {
            in_loop = true;
            table.col_names.clear();
            table.rows.clear();
            // New loop resets column collection and buffered row tokens.
            if (starts_with(line, "loop_")) {
            continue;
        }

        if (!in_loop) {
            continue;
        }

            // Collect column names matching the target prefix.
            if (starts_with(line, "_")) {
            if (!table.col_names.empty()) {
                return table;
            }
            in_loop = false;
            continue;
        }

        if (starts_with(line, "_")) {
            if (!starts_with(line, prefix)) {
                if (!table.col_names.empty()) {
                    return table;
                }
                continue;
            }
            // Rows may span multiple lines; accumulate until the row is complete.
            pending.insert(pending.end(), tokens.begin(), tokens.end());
            if (!tokens.empty()) {
                string col = tokens[0].substr(prefix.size());
                table.col_names.push_back(col);
            }
            continue;
        }

        if (table.col_names.empty()) {
            continue;
        }

        if (starts_with(line, ";")) {
            in_multiline = true;
            multiline_value = line.substr(1) + "\n";
            continue;
        }

        vector<string> tokens = split_tokens(line);
        if (tokens.empty()) {
            continue;
        }

        pending.insert(pending.end(), tokens.begin(), tokens.end());
        if (pending.size() == table.col_names.size()) {
            table.rows.push_back(pending);
            pending.clear();
        }
    }

    return table;
}

unordered_map<string, string> parse_single_value_category(const vector<string> &lines, const string &prefix) {
    unordered_map<string, string> values;
    for (const string &raw : lines) {
        string line = raw;
        if (!line.empty() && line.back() == '\n') {
            line.pop_back();
        }
        if (!starts_with(line, prefix)) {
            continue;
        }
        vector<string> tokens = split_tokens(line);
        if (tokens.size() < 2) {
            continue;
        }
        string key = tokens[0].substr(prefix.size());
        values[key] = tokens[1];
    }
    return values;
}

bool is_missing(const string &value) {
    return value == "." || value == "?";
}

array<char, 5> to_char_array_5(const string &value) {
    array<char, 5> result{};
    size_t count = 0;
    for (char ch : value) {
        if (count >= 4) {
            break;
        }
        result[count++] = ch;
    }
    result[count] = '\0';
    return result;
}

array<char, 2> to_char_array_2(const string &value) {
    array<char, 2> result{};
    if (!value.empty()) {
        result[0] = value[0];
    } else {
        result[0] = ' ';
    }
    result[1] = '\0';
    return result;
}

string sanitize_value(const string &value) {
    if (is_missing(value)) {
        return "";
    }
    return trim_quotes(value);
}

string trim_whitespace(const string &value) {
    size_t start = 0;
    while (start < value.size() && isspace(static_cast<unsigned char>(value[start]))) {
        ++start;
    }
    if (start == value.size()) {
        return "";
    }
    size_t end = value.size();
    while (end > start && isspace(static_cast<unsigned char>(value[end - 1]))) {
        --end;
    }
    return value.substr(start, end - start);
}

} // namespace

Coor MMCIF_parse(const string &filename) {
    Coor coor;
    Model model;
    coor.clear();
    model.clear();

    ifstream file(filename);
    if (!file) {
        cerr << "Error: cannot open file " << filename << endl;
        return coor;
    }

    vector<string> lines;
    string line;
    while (getline(file, line)) {
        lines.push_back(line + "\n");
    }

    LoopTable atom_site = parse_loop_for_prefix(lines, "_atom_site.");
    if (atom_site.col_names.empty()) {
        return coor;
    }

    unordered_map<string, size_t> col_index;
    for (size_t i = 0; i < atom_site.col_names.size(); ++i) {
        col_index[atom_site.col_names[i]] = i;
    }

    auto get_value = [&](const vector<string> &row, const string &key) -> string {
        auto it = col_index.find(key);
        if (it == col_index.end()) {
            return "";
        }
        if (it->second >= row.size()) {
            return "";
        }
        return row[it->second];
    };

    int uniq_resid = -1;
    int old_resid = -99999999;
    string old_chain;
    int old_model = -1;

    for (const auto &row : atom_site.rows) {
        string group_pdb = sanitize_value(get_value(row, "group_PDB"));
        bool field = (group_pdb == "HETATM");

        string num_str = sanitize_value(get_value(row, "id"));
        string name = sanitize_value(get_value(row, "label_atom_id"));
        string resname = sanitize_value(get_value(row, "label_comp_id"));
        string altloc = sanitize_value(get_value(row, "label_alt_id"));
        string chain_label = sanitize_value(get_value(row, "label_asym_id"));
        string chain_auth = sanitize_value(get_value(row, "auth_asym_id"));
        string resid_str = sanitize_value(get_value(row, "auth_seq_id"));
        string insert_res = sanitize_value(get_value(row, "pdbx_PDB_ins_code"));
        string elem = sanitize_value(get_value(row, "type_symbol"));

        string model_str = sanitize_value(get_value(row, "pdbx_PDB_model_num"));

        int num = num_str.empty() ? 0 : stoi(num_str);
        int resid = resid_str.empty() ? 0 : stoi(resid_str);
        int model_num = model_str.empty() ? 1 : stoi(model_str);

        if (chain_auth.empty()) {
            chain_auth = chain_label;
        }

        // Split models based on pdbx_PDB_model_num.
        if (model_num != old_model) {
            if (model.size() > 0) {
                coor.add_Model(model);
                model.clear();
            }
            uniq_resid = -1;
            old_resid = -99999999;
            old_chain.clear();
            old_model = model_num;
        }

        // Increment unique residue when residue id or chain changes.
        if (resid != old_resid || chain_auth != old_chain) {
            ++uniq_resid;
            old_resid = resid;
            old_chain = chain_auth;
        }

        string x_str = sanitize_value(get_value(row, "Cartn_x"));
        string y_str = sanitize_value(get_value(row, "Cartn_y"));
        string z_str = sanitize_value(get_value(row, "Cartn_z"));
        float x = x_str.empty() ? 0.0f : stof(x_str);
        float y = y_str.empty() ? 0.0f : stof(y_str);
        float z = z_str.empty() ? 0.0f : stof(z_str);

        string occ_str = sanitize_value(get_value(row, "occupancy"));
        string beta_str = sanitize_value(get_value(row, "B_iso_or_equiv"));
        float occ = occ_str.empty() ? 0.0f : stof(occ_str);
        float beta = beta_str.empty() ? 0.0f : stof(beta_str);

        array<char, 5> name_array = to_char_array_5(name);
        array<char, 5> resname_array = to_char_array_5(resname);
        array<char, 5> elem_array = to_char_array_5(elem);
        array<char, 2> alterloc_array = to_char_array_2(altloc);
        array<char, 2> chain_array = to_char_array_2(chain_label);
        array<char, 2> insert_array = to_char_array_2(insert_res);

        model.addAtom(
            num,
            name_array,
            resname_array,
            resid,
            chain_array,
            x, y, z, occ, beta,
            alterloc_array,
            elem_array,
            insert_array,
            field,
            uniq_resid
        );
    }

    if (model.size() > 0) {
        coor.add_Model(model);
    }

    auto cell_values = parse_single_value_category(lines, "_cell.");
    if (!cell_values.empty() &&
        cell_values.find("length_a") != cell_values.end() &&
        cell_values.find("length_b") != cell_values.end() &&
        cell_values.find("length_c") != cell_values.end() &&
        cell_values.find("angle_alpha") != cell_values.end() &&
        cell_values.find("angle_beta") != cell_values.end() &&
        cell_values.find("angle_gamma") != cell_values.end()) {
        float a = stof(cell_values["length_a"]);
        float b = stof(cell_values["length_b"]);
        float c = stof(cell_values["length_c"]);
        float alpha = stof(cell_values["angle_alpha"]);
        float beta = stof(cell_values["angle_beta"]);
        float gamma = stof(cell_values["angle_gamma"]);
        int z = 1;
        auto z_it = cell_values.find("Z_PDB");
        if (z_it != cell_values.end() && !z_it->second.empty()) {
            z = stoi(z_it->second);
        }

        string sgroup = "P 1";
        auto sym_values = parse_single_value_category(lines, "_symmetry.");
        auto sg_it = sym_values.find("space_group_name_H-M");
        if (sg_it != sym_values.end()) {
            sgroup = trim_quotes(sg_it->second);
        }

        ostringstream cryst;
        cryst << "CRYST1"
              << setw(9) << fixed << setprecision(3) << a
              << setw(9) << fixed << setprecision(3) << b
              << setw(9) << fixed << setprecision(3) << c
              << setw(7) << fixed << setprecision(2) << alpha
              << setw(7) << fixed << setprecision(2) << beta
              << setw(7) << fixed << setprecision(2) << gamma
              << " " << left << setw(9) << sgroup
              << right << setw(3) << z << "\n";
        coor.crystal_pack.set_CRYST1_pdb(cryst.str());
    }

    LoopTable assembly_gen = parse_loop_for_prefix(lines, "_pdbx_struct_assembly_gen.");
    LoopTable oper_list = parse_loop_for_prefix(lines, "_pdbx_struct_oper_list.");

    if (!assembly_gen.col_names.empty() && !oper_list.col_names.empty()) {
        unordered_map<string, size_t> oper_index;
        for (size_t i = 0; i < oper_list.col_names.size(); ++i) {
            oper_index[oper_list.col_names[i]] = i;
        }

        unordered_map<string, vector<vector<float>>> matrices;
        auto get_oper_val = [&](const vector<string> &row, const string &key) -> string {
            auto it = oper_index.find(key);
            if (it == oper_index.end()) {
                return "";
            }
            if (it->second >= row.size()) {
                return "";
            }
            return row[it->second];
        };

        for (const auto &row : oper_list.rows) {
            string id = sanitize_value(get_oper_val(row, "id"));
            if (id.empty()) {
                continue;
            }
            vector<vector<float>> matrix_rows;
            for (int i = 1; i <= 3; ++i) {
                vector<float> row_vals;
                for (int j = 1; j <= 3; ++j) {
                    string key = "matrix[" + to_string(i) + "][" + to_string(j) + "]";
                    row_vals.push_back(stof(sanitize_value(get_oper_val(row, key))));
                }
                string vec_key = "vector[" + to_string(i) + "]";
                row_vals.push_back(stof(sanitize_value(get_oper_val(row, vec_key))));
                matrix_rows.push_back(row_vals);
            }
            matrices[id] = matrix_rows;
        }

        unordered_map<string, size_t> gen_index;
        for (size_t i = 0; i < assembly_gen.col_names.size(); ++i) {
            gen_index[assembly_gen.col_names[i]] = i;
        }

        auto get_gen_val = [&](const vector<string> &row, const string &key) -> string {
            auto it = gen_index.find(key);
            if (it == gen_index.end()) {
                return "";
            }
            if (it->second >= row.size()) {
                return "";
            }
            return row[it->second];
        };

        if (!assembly_gen.rows.empty()) {
            const auto &row = assembly_gen.rows.front();
            string chains_str = sanitize_value(get_gen_val(row, "asym_id_list"));
            string oper_expr = sanitize_value(get_gen_val(row, "oper_expression"));

            vector<string> chains;
            if (!chains_str.empty()) {
                string token;
                istringstream ss(chains_str);
                while (getline(ss, token, ',')) {
                    token.erase(token.begin(), find_if(token.begin(), token.end(), [](unsigned char ch) { return !isspace(ch); }));
                    token.erase(find_if(token.rbegin(), token.rend(), [](unsigned char ch) { return !isspace(ch); }).base(), token.end());
                    if (!token.empty()) {
                        chains.push_back(token);
                    }
                }
            }

            vector<string> oper_ids;
            if (!oper_expr.empty()) {
                string token;
                istringstream ss(oper_expr);
                while (getline(ss, token, ',')) {
                    token.erase(remove(token.begin(), token.end(), '('), token.end());
                    token.erase(remove(token.begin(), token.end(), ')'), token.end());
                    token.erase(token.begin(), find_if(token.begin(), token.end(), [](unsigned char ch) { return !isspace(ch); }));
                    token.erase(find_if(token.rbegin(), token.rend(), [](unsigned char ch) { return !isspace(ch); }).base(), token.end());
                    if (!token.empty()) {
                        oper_ids.push_back(token);
                    }
                }
            }

            coor.transformation.clear();
            coor.transformation.set_chains(chains);
            for (const auto &oper_id : oper_ids) {
                auto it = matrices.find(oper_id);
                if (it == matrices.end()) {
                    continue;
                }
                for (const auto &row_vals : it->second) {
                    coor.transformation.add_matrix_row(row_vals);
                }
            }
        }
    }

    return coor;
}

string get_mmcif_string(const Coor &coor) {
    ostringstream oss;
    oss << "data_XXXX\n#\n";
    // Emit a minimal atom_site loop for round-trip support.
    oss << "loop_\n"
        << "_atom_site.group_PDB\n"
        << "_atom_site.id\n"
        << "_atom_site.type_symbol\n"
        << "_atom_site.label_atom_id\n"
        << "_atom_site.label_alt_id\n"
        << "_atom_site.label_comp_id\n"
        << "_atom_site.label_asym_id\n"
        << "_atom_site.label_entity_id\n"
        << "_atom_site.label_seq_id\n"
        << "_atom_site.pdbx_PDB_ins_code\n"
        << "_atom_site.Cartn_x\n"
        << "_atom_site.Cartn_y\n"
        << "_atom_site.Cartn_z\n"
        << "_atom_site.occupancy\n"
        << "_atom_site.B_iso_or_equiv\n"
        << "_atom_site.pdbx_formal_charge\n"
        << "_atom_site.auth_seq_id\n"
        << "_atom_site.auth_comp_id\n"
        << "_atom_site.auth_asym_id\n"
        << "_atom_site.auth_atom_id\n"
        << "_atom_site.pdbx_PDB_model_num\n";

    for (size_t model_index = 0; model_index < coor.model_size(); ++model_index) {
        const Model &model = coor.get_Models(model_index);
        for (size_t i = 0; i < model.size(); ++i) {
            string group = model.get_field()[i] ? "HETATM" : "ATOM";
            string name = trim_whitespace(string(model.get_name()[i].data()));
            string resname = trim_whitespace(string(model.get_resname()[i].data()));
            string chain = trim_whitespace(string(model.get_chain()[i].data()));
            string altloc = trim_whitespace(string(model.get_alterloc()[i].data()));
            string insertres = trim_whitespace(string(model.get_insertres()[i].data()));
            string elem = trim_whitespace(string(model.get_elem()[i].data()));

            if (altloc.empty()) {
                altloc = ".";
            }
            if (insertres.empty()) {
                insertres = "?";
            }
            if (elem.empty()) {
                elem = "?";
            }

            oss << group << " "
                << model.get_num()[i] << " "
                << elem << " "
                << name << " "
                << altloc << " "
                << resname << " "
                << chain << " "
                << "1" << " "
                << model.get_uniqresid()[i] + 1 << " "
                << insertres << " "
                << fixed << setprecision(3) << model.get_x()[i] << " "
                << fixed << setprecision(3) << model.get_y()[i] << " "
                << fixed << setprecision(3) << model.get_z()[i] << " "
                << fixed << setprecision(2) << model.get_occ()[i] << " "
                << fixed << setprecision(2) << model.get_beta()[i] << " "
                << "?" << " "
                << model.get_resid()[i] << " "
                << resname << " "
                << chain << " "
                << name << " "
                << (model_index + 1) << "\n";
        }
    }

    oss << "#\n";
    return oss.str();
}

bool MMCIF_write(const Coor &coor, const string &filename) {
    ofstream file(filename);
    if (!file) {
        cerr << "Error: cannot open file " << filename << endl;
        return false;
    }
    file << get_mmcif_string(coor);
    return true;
}
