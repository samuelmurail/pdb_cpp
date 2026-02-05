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

#include <stdexcept>
#include <sstream>
#include <algorithm>
#include <array>
#include <cmath>

using std::string;
using std::vector;

namespace {

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
    const std::vector<std::string> &chain_2)
{
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
