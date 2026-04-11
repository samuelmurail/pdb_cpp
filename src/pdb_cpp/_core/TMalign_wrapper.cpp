// Interface types for TM-align integration
#include "TMalign_iface.h"

// Use the vendored USalign Kabsch implementation in this translation unit.
#include "usalign/Kabsch.h"

// Disable the local Kabsch implementation in geom.h for this
// translation unit to avoid conflicting with the implementation
// provided by USalign's Kabsch.h.
#define PDBCPP_DISABLE_LOCAL_KABSCH

#include "Coor.h"
#undef PDBCPP_DISABLE_LOCAL_KABSCH


#include "Model.h"

// Use the vendored USalign implementation of TM-align
#include "usalign/TMalign.h"
// Multimer alignment helpers from USalign
#include "usalign/MMalign.h"

#include <stdexcept>
#include <cstdio>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <utility>
#include <cstdlib>
#include <atomic>

using std::string;
using std::vector;

namespace {

string array_to_string(const std::array<char, 2> &value)
{
    string out;
    for (char c : value)
    {
        if (c != '\0' && c != ' ') out.push_back(c);
    }
    return out;
}

string array_to_string(const std::array<char, 5> &value)
{
    string out;
    for (char c : value)
    {
        if (c != '\0') out.push_back(c);
    }
    return out;
}

string trim_string(string value)
{
    const size_t first = value.find_first_not_of(' ');
    if (first == string::npos) return "";
    const size_t last = value.find_last_not_of(' ');
    return value.substr(first, last - first + 1);
}

string padded_resname(const std::array<char, 5> &value)
{
    string out = trim_string(array_to_string(value));
    if (out.size() < 3) out.insert(out.begin(), 3 - out.size(), ' ');
    if (out.size() > 3) out.resize(3);
    return out;
}

string atom_name_for_pdb(const std::array<char, 5> &value)
{
    const string trimmed = trim_string(array_to_string(value));
    if (trimmed == "CA") return " CA ";
    if (trimmed == "C3'") return " C3'";

    string out = trimmed;
    if (out.size() < 4) out.insert(out.begin(), 4 - out.size(), ' ');
    if (out.size() > 4) out.resize(4);
    return out;
}

bool chain_is_requested(const string &chain_id, const vector<string> &chain2parse)
{
    if (chain2parse.empty()) return true;
    if (chain_id.empty())
    {
        return std::find(chain2parse.begin(), chain2parse.end(), "_") != chain2parse.end();
    }
    return std::find(chain2parse.begin(), chain2parse.end(), chain_id) != chain2parse.end();
}

bool select_alignment_atom(const string &atom_name)
{
    return atom_name == " CA " || atom_name == " C3'";
}

vector<vector<string>> build_pdb_lines_from_coor(
    const Coor &coor,
    const vector<string> &chain2parse,
    vector<int> &mol_vec)
{
    vector<vector<string>> pdb_lines;

    if (coor.model_size() == 0) return pdb_lines;

    const Model model = coor.get_Models(static_cast<int>(coor.get_active_model()));
    const auto &names = model.get_name();
    const auto &resnames = model.get_resname();
    const auto &chains = model.get_chain();
    const auto &resids = model.get_resid();
    const auto &insertres = model.get_insertres();
    const auto &xs = model.get_x();
    const auto &ys = model.get_y();
    const auto &zs = model.get_z();
    const auto &fields = model.get_field();
    const auto &alterlocs = model.get_alterloc();

    string previous_chain;
    int serial = 0;

    for (size_t index = 0; index < model.size(); ++index)
    {
        if (fields[index]) continue; // skip HETATM, matching het_opt=0

        const string chain_id = array_to_string(chains[index]);
        if (!chain_is_requested(chain_id, chain2parse)) continue;

        const char altloc = alterlocs[index][0];
        if (!(altloc == '\0' || altloc == ' ' || altloc == 'A')) continue;

        const string atom_name = atom_name_for_pdb(names[index]);
        const string resname = padded_resname(resnames[index]);
        if (!select_alignment_atom(atom_name)) continue;

        if (pdb_lines.empty() || previous_chain != chain_id)
        {
            pdb_lines.emplace_back();
            mol_vec.push_back(0);
            previous_chain = chain_id;
        }

        mol_vec.back() += (atom_name == " C3'") ? 1 : -1;

        std::ostringstream resi_stream;
        resi_stream << resids[index];
        const char insertion_code = insertres[index][0];
        if (insertion_code != '\0' && insertion_code != ' ')
        {
            resi_stream << insertion_code;
        }
        else
        {
            resi_stream << ' ';
        }
        string resi = resi_stream.str();
        if (resi.size() > 5) resi = resi.substr(0, 5);

        std::ostringstream line;
        line << "ATOM  "
             << std::setw(5) << (++serial) << " "
             << atom_name << " "
             << resname << " "
             << std::setw(1) << (chain_id.empty() ? ' ' : chain_id[0])
             << std::setw(5) << resi
             << "   "
             << std::setw(8) << std::fixed << std::setprecision(3) << xs[index]
             << std::setw(8) << std::fixed << std::setprecision(3) << ys[index]
             << std::setw(8) << std::fixed << std::setprecision(3) << zs[index];
        pdb_lines.back().push_back(line.str());
    }

    return pdb_lines;
}

void parse_coor_chains(
    const Coor &coor,
    vector<vector<vector<double>>> &a_vec,
    vector<vector<char>> &seq_vec,
    vector<vector<char>> &sec_vec,
    vector<int> &mol_vec,
    vector<int> &len_vec,
    int &len_aa,
    int &len_na,
    vector<string> &resi_vec,
    const vector<string> &chain2parse)
{
    vector<int> parsed_mol_vec;
    vector<int> accepted_mol_vec;
    vector<vector<string>> pdb_lines = build_pdb_lines_from_coor(coor, chain2parse, parsed_mol_vec);

    vector<double> tmp_atom_array(3, 0.0);
    const int read_resi = 2;

    for (size_t chain_i = 0; chain_i < pdb_lines.size(); ++chain_i)
    {
        int len = static_cast<int>(pdb_lines[chain_i].size());
        if (len == 0) continue;
        if (len < 3)
        {
            std::cerr << "Sequence is too short <3!" << std::endl;
            continue;
        }

        double **xa;
        NewArray(&xa, len, 3);
        char *seq = new char[len + 1];
        char *sec = new char[len + 1];
        len = read_PDB(pdb_lines[chain_i], xa, seq, resi_vec, read_resi);
        if (parsed_mol_vec[chain_i] > 0)
        {
            make_sec(seq, xa, len, sec, "auto");
        }
        else
        {
            make_sec(xa, len, sec);
        }

        vector<vector<double>> tmp_chain_array(len, tmp_atom_array);
        vector<char> tmp_seq_array(len + 1, 0);
        vector<char> tmp_sec_array(len + 1, 0);
        for (int r = 0; r < len; ++r)
        {
            tmp_chain_array[r][0] = xa[r][0];
            tmp_chain_array[r][1] = xa[r][1];
            tmp_chain_array[r][2] = xa[r][2];
            tmp_seq_array[r] = seq[r];
            tmp_sec_array[r] = sec[r];
        }

        a_vec.push_back(tmp_chain_array);
        seq_vec.push_back(tmp_seq_array);
        sec_vec.push_back(tmp_sec_array);
        len_vec.push_back(len);
        accepted_mol_vec.push_back(parsed_mol_vec[chain_i]);

        DeleteArray(&xa, len);
        delete[] seq;
        delete[] sec;
    }

    mol_vec.swap(accepted_mol_vec);
    len_aa = 0;
    len_na = 0;
    for (size_t i = 0; i < a_vec.size(); ++i)
    {
        if (mol_vec[i] > 0) len_na += len_vec[i];
        else len_aa += len_vec[i];
    }
}

// Join a list of strings with spaces (e.g. ["A","B"] -> "A B")
string join_strings(const vector<string> &values)
{
    if (values.empty()) return "";
    std::ostringstream oss;
    for (size_t i = 0; i < values.size(); ++i)
    {
        if (i > 0) oss << ' ';
        oss << values[i];
    }
    return oss.str();
}

// Build a selection string for backbone CA atoms of the requested chains
string build_ca_selection(const vector<string> &chains)
{
    const string chain_str = join_strings(chains);
    // Keep this consistent with other selections in the project
    return "chain " + chain_str +
           " and protein and name CA and not altloc B C D E F";
}

struct TempPdbFile {
    string path;

    TempPdbFile() = default;
    explicit TempPdbFile(string p) : path(std::move(p)) {}

    TempPdbFile(const TempPdbFile&) = delete;
    TempPdbFile& operator=(const TempPdbFile&) = delete;

    TempPdbFile(TempPdbFile&& other) noexcept : path(std::move(other.path)) {
        other.path.clear();
    }

    TempPdbFile& operator=(TempPdbFile&& other) noexcept {
        if (this != &other) {
            if (!path.empty()) {
                std::remove(path.c_str());
            }
            path = std::move(other.path);
            other.path.clear();
        }
        return *this;
    }

    ~TempPdbFile() {
        if (!path.empty()) {
            std::remove(path.c_str());
        }
    }
};

string make_temp_pdb_path()
{
    // Use an atomic sequence number for uniqueness — avoids platform-specific
    // mkstemp (POSIX) / _mktemp_s (MSVC) which behave differently on Windows CI.
    static std::atomic<unsigned long long> seq{0};
    unsigned long long n = seq.fetch_add(1, std::memory_order_relaxed);
#ifdef _WIN32
    const char *tmpdir = std::getenv("TEMP");
    if (!tmpdir || !tmpdir[0]) tmpdir = std::getenv("TMP");
    if (!tmpdir || !tmpdir[0]) tmpdir = "C:\\Temp";
    return string(tmpdir) + "\\pdbcpp_tmp_" + std::to_string(n) + ".pdb";
#else
    const char *tmpdir = std::getenv("TMPDIR");
    if (!tmpdir || !tmpdir[0]) tmpdir = "/tmp";
    return string(tmpdir) + "/pdbcpp_tmp_" + std::to_string(n) + ".pdb";
#endif
}

TempPdbFile write_temp_pdb(const Coor &coor)
{
    TempPdbFile temp{make_temp_pdb_path()};
    if (!coor.write(temp.path)) {
        throw std::runtime_error("Failed to write temporary PDB file: " + temp.path);
    }
    return temp;
}

TMalignResult mmalign_collect_result(
    const vector<vector<vector<double> > > &xa_vec,
    const vector<vector<vector<double> > > &ya_vec,
    const vector<vector<char> > &seqx_vec,
    const vector<vector<char> > &seqy_vec,
    const vector<vector<char> > &secx_vec,
    const vector<vector<char> > &secy_vec,
    const vector<int> &mol_vec1,
    const vector<int> &mol_vec2,
    const vector<int> &xlen_vec,
    const vector<int> &ylen_vec,
    int len_aa,
    int len_na,
    int chain1_num,
    int chain2_num,
    vector<vector<string> > &seqxA_mat,
    vector<vector<string> > &seqyA_mat,
    int *assign1_list,
    int *assign2_list,
    vector<string> &sequence,
    double d0_scale,
    bool fast_opt)
{
    int xlen = 0;
    int ylen = 0;
    for (int i = 0; i < chain1_num; ++i) xlen += xlen_vec[i];
    for (int j = 0; j < chain2_num; ++j) ylen += ylen_vec[j];

    if (xlen <= 3 || ylen <= 3) {
        throw std::runtime_error("MMalign requires at least 3 CA atoms per complex");
    }

    char *seqx = new char[xlen + 1];
    char *seqy = new char[ylen + 1];
    char *secx = new char[xlen + 1];
    char *secy = new char[ylen + 1];
    double **xa;
    double **ya;
    NewArray(&xa, xlen, 3);
    NewArray(&ya, ylen, 3);

    int mol_type = copy_chain_pair_data(
        xa_vec, ya_vec, seqx_vec, seqy_vec, secx_vec, secy_vec,
        mol_vec1, mol_vec2, xlen_vec, ylen_vec,
        xa, ya, seqx, seqy, secx, secy,
        chain1_num, chain2_num,
        seqxA_mat, seqyA_mat, assign1_list, assign2_list, sequence);

    double t0[3] = {0.0, 0.0, 0.0};
    double u0[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    double TM1 = 0.0, TM2 = 0.0, TM3 = 0.0, TM4 = 0.0, TM5 = 0.0;
    double d0_0 = 0.0, TM_0 = 0.0;
    double d0A = 0.0, d0B = 0.0, d0u = 0.0, d0a = 0.0, d0_out = 5.0;
    string seqM;
    string seqxA;
    string seqyA;
    vector<double> do_vec;
    double rmsd0 = 0.0;
    int L_ali = 0;
    double Liden = 0.0;
    double TM_ali = 0.0;
    double rmsd_ali = 0.0;
    int n_ali = 0;
    int n_ali8 = 0;
    double Lnorm_ass = len_aa + len_na;

    TMalign_main(
        xa, ya, seqx, seqy, secx, secy,
        t0, u0, TM1, TM2, TM3, TM4, TM5, d0_0, TM_0,
        d0A, d0B, d0u, d0a, d0_out,
        seqM, seqxA, seqyA, do_vec,
        rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
        xlen, ylen, sequence, Lnorm_ass, d0_scale,
        3, 0, false, false, fast_opt, mol_type, -1);

    TMalignResult result;
    result.TM1 = TM1;
    result.TM2 = TM2;
    result.TM_ali = TM_0;
    result.rmsd = rmsd0;
    result.L_ali = n_ali8;
    result.Liden = Liden;
    result.seqM = seqM;
    result.seqxA = seqxA;
    result.seqyA = seqyA;

    delete[] seqx;
    delete[] seqy;
    delete[] secx;
    delete[] secy;
    DeleteArray(&xa, xlen);
    DeleteArray(&ya, ylen);
    do_vec.clear();

    return result;
}

TMalignResult mmalign_from_files(
    const string &xname,
    const string &yname,
    const vector<string> &chain2parse1,
    const vector<string> &chain2parse2)
{
    vector<vector<vector<double> > > xa_vec;
    vector<vector<vector<double> > > ya_vec;
    vector<vector<char> > seqx_vec;
    vector<vector<char> > seqy_vec;
    vector<vector<char> > secx_vec;
    vector<vector<char> > secy_vec;
    vector<int> mol_vec1;
    vector<int> mol_vec2;
    vector<string> chainID_list1;
    vector<string> chainID_list2;
    vector<int> xlen_vec;
    vector<int> ylen_vec;
    int xlen_aa = 0;
    int xlen_na = 0;
    int ylen_aa = 0;
    int ylen_na = 0;
    vector<string> resi_vec1;
    vector<string> resi_vec2;

    vector<string> chain1_list = {xname};
    vector<string> chain2_list = {yname};
    vector<string> model2parse1;
    vector<string> model2parse2;
    vector<string> sequence;

    const double d0_scale = 0.0;
    const bool a_opt = false;
    const bool d_opt = false;
    const double TMcut = -1.0;
    const int infmt1_opt = -1;
    const int infmt2_opt = -1;
    const int ter_opt = 1;
    const int split_opt = 2;
    bool fast_opt = false;
    const int mirror_opt = 0;
    const int het_opt = 0;
    const bool autojustify = true;
    const string atom_opt = "auto";
    const string mol_opt = "auto";
    const int byresi_opt = 0;
    const bool se_opt = false;

    parse_chain_list(chain1_list, xa_vec, seqx_vec, secx_vec, mol_vec1,
        xlen_vec, chainID_list1, ter_opt, split_opt, mol_opt, infmt1_opt,
        atom_opt, autojustify, mirror_opt, het_opt, xlen_aa, xlen_na, 0,
        resi_vec1, chain2parse1, model2parse1);
    if (xa_vec.empty()) {
        throw std::runtime_error("No chains parsed from complex 1");
    }
    parse_chain_list(chain2_list, ya_vec, seqy_vec, secy_vec, mol_vec2,
        ylen_vec, chainID_list2, ter_opt, split_opt, mol_opt, infmt2_opt,
        atom_opt, autojustify, 0, het_opt, ylen_aa, ylen_na, 0,
        resi_vec2, chain2parse2, model2parse2);
    if (ya_vec.empty()) {
        throw std::runtime_error("No chains parsed from complex 2");
    }

    int len_aa = getmin(xlen_aa, ylen_aa);
    int len_na = getmin(xlen_na, ylen_na);
    if (a_opt) {
        len_aa = (xlen_aa + ylen_aa) / 2;
        len_na = (xlen_na + ylen_na) / 2;
    }

    const int chain1_num = xa_vec.size();
    const int chain2_num = ya_vec.size();

    if (chain1_num == 1 && chain2_num == 1) {
        int xlen = xlen_vec[0];
        int ylen = ylen_vec[0];
        char *seqx = new char[xlen + 1];
        char *seqy = new char[ylen + 1];
        char *secx = new char[xlen + 1];
        char *secy = new char[ylen + 1];
        double **xa;
        double **ya;
        NewArray(&xa, xlen, 3);
        NewArray(&ya, ylen, 3);
        copy_chain_data(xa_vec[0], seqx_vec[0], secx_vec[0], xlen, xa, seqx, secx);
        copy_chain_data(ya_vec[0], seqy_vec[0], secy_vec[0], ylen, ya, seqy, secy);

        double t0[3] = {0.0, 0.0, 0.0};
        double u0[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
        double TM1 = 0.0, TM2 = 0.0, TM3 = 0.0, TM4 = 0.0, TM5 = 0.0;
        double d0_0 = 0.0, TM_0 = 0.0;
        double d0A = 0.0, d0B = 0.0, d0u = 0.0, d0a = 0.0, d0_out = 5.0;
        string seqM;
        string seqxA;
        string seqyA;
        vector<double> do_vec;
        double rmsd0 = 0.0;
        int L_ali = 0;
        double Liden = 0.0;
        double TM_ali = 0.0;
        double rmsd_ali = 0.0;
        int n_ali = 0;
        int n_ali8 = 0;

        TMalign_main(xa, ya, seqx, seqy, secx, secy,
            t0, u0, TM1, TM2, TM3, TM4, TM5,
            d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
            seqM, seqxA, seqyA, do_vec,
            rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
            xlen, ylen, sequence, 0, d0_scale,
            0, a_opt, false, d_opt, fast_opt,
            mol_vec1[0] + mol_vec2[0], TMcut);

        TMalignResult result;
        result.TM1 = TM1;
        result.TM2 = TM2;
        result.TM_ali = TM_0;
        result.rmsd = rmsd0;
        result.L_ali = n_ali8;
        result.Liden = Liden;
        result.seqM = seqM;
        result.seqxA = seqxA;
        result.seqyA = seqyA;

        delete[] seqx;
        delete[] seqy;
        delete[] secx;
        delete[] secy;
        DeleteArray(&xa, xlen);
        DeleteArray(&ya, ylen);
        do_vec.clear();

        return result;
    }

    vector<string> tmp_str_vec(chain2_num, "");
    double **TMave_mat;
    double **ut_mat;
    NewArray(&TMave_mat, chain1_num, chain2_num);
    NewArray(&ut_mat, chain1_num * chain2_num, 4 * 3);
    vector<vector<string> > seqxA_mat(chain1_num, tmp_str_vec);
    vector<vector<string> > seqM_mat(chain1_num, tmp_str_vec);
    vector<vector<string> > seqyA_mat(chain1_num, tmp_str_vec);

    double maxTMmono = -1.0;
    int maxTMmono_i = 0;
    int maxTMmono_j = 0;

    if (len_aa + len_na > 500) fast_opt = true;

    for (int i = 0; i < chain1_num; ++i)
    {
        int xlen = xlen_vec[i];
        if (xlen < 3)
        {
            for (int j = 0; j < chain2_num; ++j) TMave_mat[i][j] = -1;
            continue;
        }

        char *seqx = new char[xlen + 1];
        char *secx = new char[xlen + 1];
        double **xa;
        NewArray(&xa, xlen, 3);
        copy_chain_data(xa_vec[i], seqx_vec[i], secx_vec[i], xlen, xa, seqx, secx);

        for (int j = 0; j < chain2_num; ++j)
        {
            int ut_idx = i * chain2_num + j;
            for (int ui = 0; ui < 4; ++ui)
                for (int uj = 0; uj < 3; ++uj)
                    ut_mat[ut_idx][ui * 3 + uj] = 0;
            ut_mat[ut_idx][0] = 1;
            ut_mat[ut_idx][4] = 1;
            ut_mat[ut_idx][8] = 1;

            if (mol_vec1[i] * mol_vec2[j] < 0)
            {
                TMave_mat[i][j] = -1;
                continue;
            }

            int ylen = ylen_vec[j];
            if (ylen < 3)
            {
                TMave_mat[i][j] = -1;
                continue;
            }

            char *seqy = new char[ylen + 1];
            char *secy = new char[ylen + 1];
            double **ya;
            NewArray(&ya, ylen, 3);
            copy_chain_data(ya_vec[j], seqy_vec[j], secy_vec[j], ylen, ya, seqy, secy);

            double t0[3], u0[3][3];
            double TM1 = 0.0, TM2 = 0.0, TM3 = 0.0, TM4 = 0.0, TM5 = 0.0;
            double d0_0 = 0.0, TM_0 = 0.0;
            double d0A = 0.0, d0B = 0.0, d0u = 0.0, d0a = 0.0;
            double d0_out = 5.0;
            string seqM;
            string seqxA;
            string seqyA;
            double rmsd0 = 0.0;
            int L_ali = 0;
            double Liden = 0.0;
            double TM_ali = 0.0;
            double rmsd_ali = 0.0;
            int n_ali = 0;
            int n_ali8 = 0;
            vector<double> do_vec;

            int Lnorm_tmp = len_aa;
            if (mol_vec1[i] + mol_vec2[j] > 0) Lnorm_tmp = len_na;

            TMalign_main(xa, ya, seqx, seqy, secx, secy,
                t0, u0, TM1, TM2, TM3, TM4, TM5,
                d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                seqM, seqxA, seqyA, do_vec,
                rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                xlen, ylen, sequence, Lnorm_tmp, d0_scale,
                0, false, true, false, fast_opt,
                mol_vec1[i] + mol_vec2[j], TMcut);

            for (int ui = 0; ui < 3; ++ui)
                for (int uj = 0; uj < 3; ++uj)
                    ut_mat[ut_idx][ui * 3 + uj] = u0[ui][uj];
            for (int uj = 0; uj < 3; ++uj) ut_mat[ut_idx][9 + uj] = t0[uj];
            seqxA_mat[i][j] = seqxA;
            seqyA_mat[i][j] = seqyA;
            TMave_mat[i][j] = TM4 * Lnorm_tmp;
            if (TMave_mat[i][j] > maxTMmono)
            {
                maxTMmono = TMave_mat[i][j];
                maxTMmono_i = i;
                maxTMmono_j = j;
            }

            delete[] seqy;
            delete[] secy;
            DeleteArray(&ya, ylen);
            do_vec.clear();
        }

        delete[] seqx;
        delete[] secx;
        DeleteArray(&xa, xlen);
    }

    int *assign1_list = new int[chain1_num];
    int *assign2_list = new int[chain2_num];
    double total_score = enhanced_greedy_search(TMave_mat, assign1_list,
        assign2_list, chain1_num, chain2_num);
    if (total_score <= 0) {
        throw std::runtime_error("No assignable chain mapping found");
    }

    int aln_chain_num = count_assign_pair(assign1_list, chain1_num);
    bool is_oligomer = (aln_chain_num >= 3);
    if (aln_chain_num == 2 && !se_opt)
    {
        int na_chain_num1, na_chain_num2, aa_chain_num1, aa_chain_num2;
        count_na_aa_chain_num(na_chain_num1, aa_chain_num1, mol_vec1);
        count_na_aa_chain_num(na_chain_num2, aa_chain_num2, mol_vec2);

        if (na_chain_num1 == 1 && na_chain_num2 == 1 &&
            aa_chain_num1 == 1 && aa_chain_num2 == 1)
        {
            is_oligomer = false;
        }
        else if ((getmin(na_chain_num1, na_chain_num2) == 0 &&
                    aa_chain_num1 == 2 && aa_chain_num2 == 2) ||
                 (getmin(aa_chain_num1, aa_chain_num2) == 0 &&
                    na_chain_num1 == 2 && na_chain_num2 == 2))
        {
            adjust_dimer_assignment(xa_vec, ya_vec, xlen_vec, ylen_vec, mol_vec1,
                mol_vec2, assign1_list, assign2_list, seqxA_mat, seqyA_mat);
            is_oligomer = false;
        }
        else
        {
            is_oligomer = true;
        }
    }

    if ((aln_chain_num >= 3 || is_oligomer) && !se_opt)
    {
        double **xcentroids;
        double **ycentroids;
        NewArray(&xcentroids, chain1_num, 3);
        NewArray(&ycentroids, chain2_num, 3);
        double d0MM = getmin(
            calculate_centroids(xa_vec, chain1_num, xcentroids),
            calculate_centroids(ya_vec, chain2_num, ycentroids));

        homo_refined_greedy_search(TMave_mat, assign1_list,
            assign2_list, chain1_num, chain2_num, xcentroids,
            ycentroids, d0MM, len_aa + len_na, ut_mat);
        hetero_refined_greedy_search(TMave_mat, assign1_list,
            assign2_list, chain1_num, chain2_num, xcentroids,
            ycentroids, d0MM, len_aa + len_na);

        DeleteArray(&xcentroids, chain1_num);
        DeleteArray(&ycentroids, chain2_num);
    }

    int init_pair_num = count_assign_pair(assign1_list, chain1_num);
    int *assign1_init = new int[chain1_num];
    int *assign2_init = new int[chain2_num];
    double **TMave_init;
    NewArray(&TMave_init, chain1_num, chain2_num);
    vector<vector<string> > seqxA_init(chain1_num, tmp_str_vec);
    vector<vector<string> > seqyA_init(chain1_num, tmp_str_vec);
    vector<string> sequence_init;
    copy_chain_assign_data(chain1_num, chain2_num, sequence_init,
        seqxA_mat, seqyA_mat, assign1_list, assign2_list, TMave_mat,
        seqxA_init, seqyA_init, assign1_init, assign2_init, TMave_init);

    double max_total_score = 0.0;
    int max_iter = 5 - static_cast<int>((len_aa + len_na) / 200);
    if (max_iter < 2) max_iter = 2;
    if (!se_opt)
    {
        double **xa;
        double **ya;
        char *seqx;
        char *seqy;
        char *secx;
        char *secy;
        std::map<int, int> chainmap;
        MMalign_iter(max_total_score, max_iter, xa_vec, ya_vec,
            seqx_vec, seqy_vec, secx_vec, secy_vec, mol_vec1, mol_vec2, xlen_vec,
            ylen_vec, xa, ya, seqx, seqy, secx, secy, len_aa, len_na, chain1_num,
            chain2_num, TMave_mat, seqxA_mat, seqyA_mat, assign1_list, assign2_list,
            sequence, d0_scale, fast_opt, chainmap, byresi_opt);
    }

    if (byresi_opt == 0 && max_total_score < maxTMmono)
    {
        copy_chain_assign_data(chain1_num, chain2_num, sequence,
            seqxA_init, seqyA_init, assign1_init, assign2_init, TMave_init,
            seqxA_mat, seqyA_mat, assign1_list, assign2_list, TMave_mat);
        for (int i = 0; i < chain1_num; ++i)
        {
            if (i != maxTMmono_i) assign1_list[i] = -1;
            else assign1_list[i] = maxTMmono_j;
        }
        for (int j = 0; j < chain2_num; ++j)
        {
            if (j != maxTMmono_j) assign2_list[j] = -1;
            else assign2_list[j] = maxTMmono_i;
        }
        sequence.clear();
        sequence.push_back(seqxA_mat[maxTMmono_i][maxTMmono_j]);
        sequence.push_back(seqyA_mat[maxTMmono_i][maxTMmono_j]);
    }

    double max_total_score_cross = max_total_score;
    if (byresi_opt == 0 && len_aa + len_na < 10000)
    {
        double **xa;
        double **ya;
        char *seqx;
        char *seqy;
        char *secx;
        char *secy;
        MMalign_dimer(max_total_score_cross, xa_vec, ya_vec, seqx_vec, seqy_vec,
            secx_vec, secy_vec, mol_vec1, mol_vec2, xlen_vec, ylen_vec,
            xa, ya, seqx, seqy, secx, secy, len_aa, len_na, chain1_num, chain2_num,
            TMave_init, seqxA_init, seqyA_init, assign1_init, assign2_init,
            sequence_init, d0_scale, fast_opt);
        if (max_total_score_cross > max_total_score)
        {
            max_total_score = max_total_score_cross;
            copy_chain_assign_data(chain1_num, chain2_num, sequence,
                seqxA_init, seqyA_init, assign1_init, assign2_init, TMave_init,
                seqxA_mat, seqyA_mat, assign1_list, assign2_list, TMave_mat);
        }
    }

    TMalignResult result = mmalign_collect_result(
        xa_vec, ya_vec, seqx_vec, seqy_vec, secx_vec, secy_vec,
        mol_vec1, mol_vec2, xlen_vec, ylen_vec, len_aa, len_na,
        chain1_num, chain2_num, seqxA_mat, seqyA_mat,
        assign1_list, assign2_list, sequence, d0_scale, fast_opt);

    delete [] assign1_list;
    delete [] assign2_list;
    DeleteArray(&TMave_mat, chain1_num);
    DeleteArray(&ut_mat, chain1_num * chain2_num);
    delete [] assign1_init;
    delete [] assign2_init;
    DeleteArray(&TMave_init, chain1_num);

    return result;
}

TMalignResult mmalign_from_coor(
    const Coor &coor_1,
    const Coor &coor_2,
    const vector<string> &chain2parse1,
    const vector<string> &chain2parse2)
{
    vector<vector<vector<double>>> xa_vec;
    vector<vector<vector<double>>> ya_vec;
    vector<vector<char>> seqx_vec;
    vector<vector<char>> seqy_vec;
    vector<vector<char>> secx_vec;
    vector<vector<char>> secy_vec;
    vector<int> mol_vec1;
    vector<int> mol_vec2;
    vector<int> xlen_vec;
    vector<int> ylen_vec;
    int xlen_aa = 0;
    int xlen_na = 0;
    int ylen_aa = 0;
    int ylen_na = 0;
    vector<string> resi_vec1;
    vector<string> resi_vec2;

    parse_coor_chains(coor_1, xa_vec, seqx_vec, secx_vec, mol_vec1,
        xlen_vec, xlen_aa, xlen_na, resi_vec1, chain2parse1);
    if (xa_vec.empty()) {
        throw std::runtime_error("No chains parsed from complex 1");
    }

    parse_coor_chains(coor_2, ya_vec, seqy_vec, secy_vec, mol_vec2,
        ylen_vec, ylen_aa, ylen_na, resi_vec2, chain2parse2);
    if (ya_vec.empty()) {
        throw std::runtime_error("No chains parsed from complex 2");
    }

    int len_aa = getmin(xlen_aa, ylen_aa);
    int len_na = getmin(xlen_na, ylen_na);
    const int chain1_num = xa_vec.size();
    const int chain2_num = ya_vec.size();

    if (chain1_num == 1 && chain2_num == 1) {
        int xlen = xlen_vec[0];
        int ylen = ylen_vec[0];
        char *seqx = new char[xlen + 1];
        char *seqy = new char[ylen + 1];
        char *secx = new char[xlen + 1];
        char *secy = new char[ylen + 1];
        double **xa;
        double **ya;
        NewArray(&xa, xlen, 3);
        NewArray(&ya, ylen, 3);
        copy_chain_data(xa_vec[0], seqx_vec[0], secx_vec[0], xlen, xa, seqx, secx);
        copy_chain_data(ya_vec[0], seqy_vec[0], secy_vec[0], ylen, ya, seqy, secy);

        double t0[3] = {0.0, 0.0, 0.0};
        double u0[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
        double TM1 = 0.0, TM2 = 0.0, TM3 = 0.0, TM4 = 0.0, TM5 = 0.0;
        double d0_0 = 0.0, TM_0 = 0.0;
        double d0A = 0.0, d0B = 0.0, d0u = 0.0, d0a = 0.0, d0_out = 5.0;
        string seqM;
        string seqxA;
        string seqyA;
        vector<double> do_vec;
        double rmsd0 = 0.0;
        int L_ali = 0;
        double Liden = 0.0;
        double TM_ali = 0.0;
        double rmsd_ali = 0.0;
        int n_ali = 0;
        int n_ali8 = 0;
        vector<string> sequence;

        TMalign_main(xa, ya, seqx, seqy, secx, secy,
            t0, u0, TM1, TM2, TM3, TM4, TM5,
            d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
            seqM, seqxA, seqyA, do_vec,
            rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
            xlen, ylen, sequence, 0, 0.0,
            0, false, true, false, false,
            mol_vec1[0] + mol_vec2[0], -1.0);

        TMalignResult result;
        result.TM1 = TM1;
        result.TM2 = TM2;
        result.TM_ali = TM_0;
        result.rmsd = rmsd0;
        result.L_ali = n_ali8;
        result.Liden = Liden;
        result.seqM = seqM;
        result.seqxA = seqxA;
        result.seqyA = seqyA;

        delete[] seqx;
        delete[] seqy;
        delete[] secx;
        delete[] secy;
        DeleteArray(&xa, xlen);
        DeleteArray(&ya, ylen);
        do_vec.clear();

        return result;
    }

    vector<string> sequence;
    vector<string> tmp_str_vec(chain2_num, "");
    double **TMave_mat;
    double **ut_mat;
    NewArray(&TMave_mat, chain1_num, chain2_num);
    NewArray(&ut_mat, chain1_num * chain2_num, 4 * 3);
    vector<vector<string>> seqxA_mat(chain1_num, tmp_str_vec);
    vector<vector<string>> seqyA_mat(chain1_num, tmp_str_vec);

    double maxTMmono = -1.0;
    int maxTMmono_i = 0;
    int maxTMmono_j = 0;
    bool fast_opt = false;
    if (len_aa + len_na > 500) fast_opt = true;

    for (int i = 0; i < chain1_num; ++i)
    {
        int xlen = xlen_vec[i];
        if (xlen < 3)
        {
            for (int j = 0; j < chain2_num; ++j) TMave_mat[i][j] = -1;
            continue;
        }

        char *seqx = new char[xlen + 1];
        char *secx = new char[xlen + 1];
        double **xa;
        NewArray(&xa, xlen, 3);
        copy_chain_data(xa_vec[i], seqx_vec[i], secx_vec[i], xlen, xa, seqx, secx);

        for (int j = 0; j < chain2_num; ++j)
        {
            int ut_idx = i * chain2_num + j;
            for (int ui = 0; ui < 4; ++ui)
                for (int uj = 0; uj < 3; ++uj)
                    ut_mat[ut_idx][ui * 3 + uj] = 0;
            ut_mat[ut_idx][0] = 1;
            ut_mat[ut_idx][4] = 1;
            ut_mat[ut_idx][8] = 1;

            if (mol_vec1[i] * mol_vec2[j] < 0)
            {
                TMave_mat[i][j] = -1;
                continue;
            }

            int ylen = ylen_vec[j];
            if (ylen < 3)
            {
                TMave_mat[i][j] = -1;
                continue;
            }

            char *seqy = new char[ylen + 1];
            char *secy = new char[ylen + 1];
            double **ya;
            NewArray(&ya, ylen, 3);
            copy_chain_data(ya_vec[j], seqy_vec[j], secy_vec[j], ylen, ya, seqy, secy);

            double t0[3], u0[3][3];
            double TM1 = 0.0, TM2 = 0.0, TM3 = 0.0, TM4 = 0.0, TM5 = 0.0;
            double d0_0 = 0.0, TM_0 = 0.0;
            double d0A = 0.0, d0B = 0.0, d0u = 0.0, d0a = 0.0;
            double d0_out = 5.0;
            string seqM;
            string seqxA;
            string seqyA;
            double rmsd0 = 0.0;
            int L_ali = 0;
            double Liden = 0.0;
            double TM_ali = 0.0;
            double rmsd_ali = 0.0;
            int n_ali = 0;
            int n_ali8 = 0;
            vector<double> do_vec;

            int Lnorm_tmp = len_aa;
            if (mol_vec1[i] + mol_vec2[j] > 0) Lnorm_tmp = len_na;

            TMalign_main(xa, ya, seqx, seqy, secx, secy,
                t0, u0, TM1, TM2, TM3, TM4, TM5,
                d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                seqM, seqxA, seqyA, do_vec,
                rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                xlen, ylen, sequence, Lnorm_tmp, 0.0,
                0, false, true, false, fast_opt,
                mol_vec1[i] + mol_vec2[j], -1.0);

            for (int ui = 0; ui < 3; ++ui)
                for (int uj = 0; uj < 3; ++uj)
                    ut_mat[ut_idx][ui * 3 + uj] = u0[ui][uj];
            for (int uj = 0; uj < 3; ++uj) ut_mat[ut_idx][9 + uj] = t0[uj];
            seqxA_mat[i][j] = seqxA;
            seqyA_mat[i][j] = seqyA;
            TMave_mat[i][j] = TM4 * Lnorm_tmp;
            if (TMave_mat[i][j] > maxTMmono)
            {
                maxTMmono = TMave_mat[i][j];
                maxTMmono_i = i;
                maxTMmono_j = j;
            }

            delete[] seqy;
            delete[] secy;
            DeleteArray(&ya, ylen);
            do_vec.clear();
        }

        delete[] seqx;
        delete[] secx;
        DeleteArray(&xa, xlen);
    }

    int *assign1_list = new int[chain1_num];
    int *assign2_list = new int[chain2_num];
    double total_score = enhanced_greedy_search(TMave_mat, assign1_list,
        assign2_list, chain1_num, chain2_num);
    if (total_score <= 0) {
        throw std::runtime_error("No assignable chain mapping found");
    }

    int aln_chain_num = count_assign_pair(assign1_list, chain1_num);
    bool is_oligomer = (aln_chain_num >= 3);
    if (aln_chain_num == 2)
    {
        int na_chain_num1, na_chain_num2, aa_chain_num1, aa_chain_num2;
        count_na_aa_chain_num(na_chain_num1, aa_chain_num1, mol_vec1);
        count_na_aa_chain_num(na_chain_num2, aa_chain_num2, mol_vec2);

        if (na_chain_num1 == 1 && na_chain_num2 == 1 &&
            aa_chain_num1 == 1 && aa_chain_num2 == 1)
        {
            is_oligomer = false;
        }
        else if ((getmin(na_chain_num1, na_chain_num2) == 0 &&
                    aa_chain_num1 == 2 && aa_chain_num2 == 2) ||
                 (getmin(aa_chain_num1, aa_chain_num2) == 0 &&
                    na_chain_num1 == 2 && na_chain_num2 == 2))
        {
            adjust_dimer_assignment(xa_vec, ya_vec, xlen_vec, ylen_vec, mol_vec1,
                mol_vec2, assign1_list, assign2_list, seqxA_mat, seqyA_mat);
            is_oligomer = false;
        }
        else
        {
            is_oligomer = true;
        }
    }

    if (aln_chain_num >= 3 || is_oligomer)
    {
        double **xcentroids;
        double **ycentroids;
        NewArray(&xcentroids, chain1_num, 3);
        NewArray(&ycentroids, chain2_num, 3);
        double d0MM = getmin(
            calculate_centroids(xa_vec, chain1_num, xcentroids),
            calculate_centroids(ya_vec, chain2_num, ycentroids));

        homo_refined_greedy_search(TMave_mat, assign1_list,
            assign2_list, chain1_num, chain2_num, xcentroids,
            ycentroids, d0MM, len_aa + len_na, ut_mat);
        hetero_refined_greedy_search(TMave_mat, assign1_list,
            assign2_list, chain1_num, chain2_num, xcentroids,
            ycentroids, d0MM, len_aa + len_na);

        DeleteArray(&xcentroids, chain1_num);
        DeleteArray(&ycentroids, chain2_num);
    }

    int *assign1_init = new int[chain1_num];
    int *assign2_init = new int[chain2_num];
    double **TMave_init;
    NewArray(&TMave_init, chain1_num, chain2_num);
    vector<string> tmp_str_vec2(chain2_num, "");
    vector<vector<string>> seqxA_init(chain1_num, tmp_str_vec2);
    vector<vector<string>> seqyA_init(chain1_num, tmp_str_vec2);
    vector<string> sequence_init;
    copy_chain_assign_data(chain1_num, chain2_num, sequence_init,
        seqxA_mat, seqyA_mat, assign1_list, assign2_list, TMave_mat,
        seqxA_init, seqyA_init, assign1_init, assign2_init, TMave_init);

    double max_total_score = 0.0;
    int max_iter = 5 - static_cast<int>((len_aa + len_na) / 200);
    if (max_iter < 2) max_iter = 2;
    {
        double **xa;
        double **ya;
        char *seqx;
        char *seqy;
        char *secx;
        char *secy;
        std::map<int, int> chainmap;
        const int byresi_opt = 0;
        MMalign_iter(max_total_score, max_iter, xa_vec, ya_vec,
            seqx_vec, seqy_vec, secx_vec, secy_vec, mol_vec1, mol_vec2, xlen_vec,
            ylen_vec, xa, ya, seqx, seqy, secx, secy, len_aa, len_na, chain1_num,
            chain2_num, TMave_mat, seqxA_mat, seqyA_mat, assign1_list, assign2_list,
            sequence, 0.0, fast_opt, chainmap, byresi_opt);
    }

    if (max_total_score < maxTMmono)
    {
        copy_chain_assign_data(chain1_num, chain2_num, sequence,
            seqxA_init, seqyA_init, assign1_init, assign2_init, TMave_init,
            seqxA_mat, seqyA_mat, assign1_list, assign2_list, TMave_mat);
        for (int i = 0; i < chain1_num; ++i)
        {
            if (i != maxTMmono_i) assign1_list[i] = -1;
            else assign1_list[i] = maxTMmono_j;
        }
        for (int j = 0; j < chain2_num; ++j)
        {
            if (j != maxTMmono_j) assign2_list[j] = -1;
            else assign2_list[j] = maxTMmono_i;
        }
        sequence.clear();
        sequence.push_back(seqxA_mat[maxTMmono_i][maxTMmono_j]);
        sequence.push_back(seqyA_mat[maxTMmono_i][maxTMmono_j]);
    }

    double max_total_score_cross = max_total_score;
    if (len_aa + len_na < 10000)
    {
        double **xa;
        double **ya;
        char *seqx;
        char *seqy;
        char *secx;
        char *secy;
        MMalign_dimer(max_total_score_cross, xa_vec, ya_vec, seqx_vec, seqy_vec,
            secx_vec, secy_vec, mol_vec1, mol_vec2, xlen_vec, ylen_vec,
            xa, ya, seqx, seqy, secx, secy, len_aa, len_na, chain1_num, chain2_num,
            TMave_init, seqxA_init, seqyA_init, assign1_init, assign2_init,
            sequence_init, 0.0, fast_opt);
        if (max_total_score_cross > max_total_score)
        {
            max_total_score = max_total_score_cross;
            copy_chain_assign_data(chain1_num, chain2_num, sequence_init,
                seqxA_init, seqyA_init, assign1_init, assign2_init, TMave_init,
                seqxA_mat, seqyA_mat, assign1_list, assign2_list, TMave_mat);
        }
    }

    TMalignResult result = mmalign_collect_result(
        xa_vec, ya_vec, seqx_vec, seqy_vec, secx_vec, secy_vec,
        mol_vec1, mol_vec2, xlen_vec, ylen_vec, len_aa, len_na,
        chain1_num, chain2_num, seqxA_mat, seqyA_mat,
        assign1_list, assign2_list, sequence, 0.0, fast_opt);

    delete[] assign1_list;
    delete[] assign2_list;
    DeleteArray(&TMave_mat, chain1_num);
    DeleteArray(&ut_mat, chain1_num * chain2_num);
    delete[] assign1_init;
    delete[] assign2_init;
    DeleteArray(&TMave_init, chain1_num);

    return result;
}

} // namespace

// secondary structure assignment for protein:
// 1->coil, 2->helix, 3->turn, 4->strand
std::vector<std::string> compute_SS(const Model &model, bool gap_in_seq)
{
    size_t j1, j2, j3, j4, j5, gap_num;
    double d13, d14, d15, d24, d25, d35;

    std::vector<int> CA_indexes = model.get_index_select("name CA");
    size_t res_num = CA_indexes.size();

    std::vector<std::array<char, 2>> uniq_chains = model.get_uniq_chain();
    std::vector<std::array<char, 2>> chain_array = model.get_chain();
    std::vector<int> resid_array = model.get_resid();

    std::array<char, 2> old_chain = chain_array[CA_indexes[2]];
    auto it = std::find(uniq_chains.begin(), uniq_chains.end(), old_chain);
    if (it == uniq_chains.end())
    {
        throw std::runtime_error("Chain not found in unique chains");
    }
    int chain_index = static_cast<int>(distance(uniq_chains.begin(), it));

    std::vector<std::string> seq_vec;
    int old_resid = resid_array[CA_indexes[2]];
    seq_vec.emplace_back("CC");

    for (size_t i = 2; i < res_num - 2; i++)
    {
        if (chain_array[CA_indexes[i]] != old_chain)
        {
            if (!seq_vec[chain_index].empty()) seq_vec[chain_index].pop_back();
            if (!seq_vec[chain_index].empty()) seq_vec[chain_index].pop_back();
            seq_vec[chain_index] += "CC";

            old_chain = chain_array[CA_indexes[i]];
            old_resid = resid_array[CA_indexes[i]];
            it = std::find(uniq_chains.begin(), uniq_chains.end(), old_chain);
            if (it == uniq_chains.end())
            {
                throw std::runtime_error("Chain not found in unique chains");
            }
            chain_index = static_cast<int>(distance(uniq_chains.begin(), it));
            seq_vec.emplace_back("CC");
            i += 2;
        }

        j1 = CA_indexes[i - 2];
        j2 = CA_indexes[i - 1];
        j3 = CA_indexes[i];
        j4 = CA_indexes[i + 1];
        j5 = CA_indexes[i + 2];
        d13 = model.distance(j1, j3);
        d14 = model.distance(j1, j4);
        d15 = model.distance(j1, j5);
        d24 = model.distance(j2, j4);
        d25 = model.distance(j2, j5);
        d35 = model.distance(j3, j5);

        if (resid_array[CA_indexes[i]] != old_resid)
        {
            if (gap_in_seq)
            {
                gap_num = resid_array[CA_indexes[i]] - old_resid;
                for (size_t j = 0; j < gap_num; ++j)
                {
                    seq_vec[chain_index] += "-";
                }
            }
            old_resid = resid_array[CA_indexes[i]];
        }

        seq_vec[chain_index] += sec_str(d13, d14, d15, d24, d25, d35);
        old_resid += 1;
    }
    seq_vec[chain_index] += "CC";

    return seq_vec;
}

std::vector<std::vector<std::string>> compute_SS(const Coor &coor, bool gap_in_seq)
{
    std::vector<std::vector<std::string>> ss_vec;

    for (size_t i = 0; i < coor.model_size(); ++i)
    {
        ss_vec.push_back(compute_SS(coor.get_Models(static_cast<int>(i)), gap_in_seq));
    }
    return ss_vec;
}

TMalignResult tmalign_CA(
    const Coor &coor_1,
    const Coor &coor_2,
    const std::vector<std::string> &chain_1,
    const std::vector<std::string> &chain_2,
    int mm)
{
    if (mm == 1) {
        return mmalign_from_coor(coor_1, coor_2, chain_1, chain_2);
    }

    // Select CA atoms for the requested chains in each structure
    const string sel_1 = build_ca_selection(chain_1);
    const string sel_2 = build_ca_selection(chain_2);

    std::vector<int> index_1 = coor_1.get_index_select(sel_1, 0);
    std::vector<int> index_2 = coor_2.get_index_select(sel_2, 0);

    const size_t xlen = index_1.size();
    const size_t ylen = index_2.size();

    if (xlen < 3 || ylen < 3)
    {
        throw std::runtime_error(
            "TMalign requires at least 3 CA atoms in each structure");
    }

    // Extract coordinates from the active model of each Coor
    const size_t model_idx_1 = coor_1.get_active_model();
    const size_t model_idx_2 = coor_2.get_active_model();

    Model model_1 = coor_1.get_Models(static_cast<int>(model_idx_1));
    Model model_2 = coor_2.get_Models(static_cast<int>(model_idx_2));

    const auto &x1 = model_1.get_x();
    const auto &y1 = model_1.get_y();
    const auto &z1 = model_1.get_z();

    const auto &x2 = model_2.get_x();
    const auto &y2 = model_2.get_y();
    const auto &z2 = model_2.get_z();

    // Allocate coordinate arrays expected by TMalign_main
    double **xa = new double *[xlen];
    double **ya = new double *[ylen];

    for (size_t i = 0; i < xlen; ++i)
    {
        xa[i] = new double[3];
        const int idx = index_1[i];
        xa[i][0] = static_cast<double>(x1[idx]);
        xa[i][1] = static_cast<double>(y1[idx]);
        xa[i][2] = static_cast<double>(z1[idx]);
    }

    for (size_t i = 0; i < ylen; ++i)
    {
        ya[i] = new double[3];
        const int idx = index_2[i];
        ya[i][0] = static_cast<double>(x2[idx]);
        ya[i][1] = static_cast<double>(y2[idx]);
        ya[i][2] = static_cast<double>(z2[idx]);
    }

    // Dummy sequence and secondary structure strings –
    // TM-align uses them mainly for heuristics; coordinates drive the fit.
    char *seqx = new char[xlen + 1];
    char *seqy = new char[ylen + 1];
    char *secx = new char[xlen + 1];
    char *secy = new char[ylen + 1];

    for (size_t i = 0; i < xlen; ++i)
    {
        seqx[i] = 'A';
        secx[i] = 'C';
    }
    for (size_t i = 0; i < ylen; ++i)
    {
        seqy[i] = 'A';
        secy[i] = 'C';
    }

    seqx[xlen] = '\0';
    seqy[ylen] = '\0';
    secx[xlen] = '\0';
    secy[ylen] = '\0';

    // Parameters and outputs for TMalign_main, following its usual defaults
    double t0[3] = {0.0, 0.0, 0.0};
    double u0[3][3] = {
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0}};

    double TM1 = 0.0, TM2 = 0.0, TM3 = 0.0, TM4 = 0.0, TM5 = 0.0;
    double d0_0 = 0.0, TM_0 = 0.0;
    double d0A = 0.0, d0B = 0.0, d0u = 0.0, d0a = 0.0, d0_out = 5.0;

    string seqM;
    string seqxA;
    string seqyA;
    vector<double> do_vec;

    double rmsd0 = 0.0;
    int L_ali = 0;
    double Liden = 0.0;
    double TM_ali = 0.0;
    double rmsd_ali = 0.0;
    int n_ali = 0;
    int n_ali8 = 0;

    vector<string> sequence; // no user-specified initial alignment
    const double Lnorm_ass = 0.0;   // use internal defaults
    const double d0_scale = 0.0;    // use internal defaults
    const int i_opt = 0;
    const int a_opt = 0;
    const bool u_opt = false;
    const bool d_opt = false;
    const bool fast_opt = false;
    const int mol_type = -1;        // protein
    const double TMcut = -1.0;      // no early cutoff

    const int status = TMalign_main(
        xa, ya,
        seqx, seqy,
        secx, secy,
        t0, u0,
        TM1, TM2, TM3, TM4, TM5,
        d0_0, TM_0,
        d0A, d0B, d0u, d0a, d0_out,
        seqM, seqxA, seqyA, do_vec,
        rmsd0, L_ali, Liden,
        TM_ali, rmsd_ali, n_ali, n_ali8,
        static_cast<int>(xlen), static_cast<int>(ylen),
        sequence, Lnorm_ass, d0_scale,
        i_opt, a_opt, u_opt, d_opt, fast_opt,
        mol_type, TMcut);

    // Clean up allocated resources
    for (size_t i = 0; i < xlen; ++i)
        delete[] xa[i];
    for (size_t i = 0; i < ylen; ++i)
        delete[] ya[i];
    delete[] xa;
    delete[] ya;

    delete[] seqx;
    delete[] seqy;
    delete[] secx;
    delete[] secy;

    if (status != 0)
    {
        throw std::runtime_error("TMalign_main did not complete a full TM-score calculation");
    }

    TMalignResult result;
    result.TM1 = TM1;
    result.TM2 = TM2;
    result.TM_ali = TM_0;
    result.rmsd = rmsd0;
    result.L_ali = n_ali8;
    result.Liden = Liden;
    result.seqM = seqM;
    result.seqxA = seqxA;
    result.seqyA = seqyA;

    return result;
}
