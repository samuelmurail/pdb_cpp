#ifndef PDBCPP_TMALIGN_IFACE_H
#define PDBCPP_TMALIGN_IFACE_H

#include <string>
#include <vector>

class Coor;
class Model;

// Result structure for TM-align based structural alignment
struct TMalignResult {
    double TM1;      // TM-score normalized by length of structure 1
    double TM2;      // TM-score normalized by length of structure 2
    double TM_ali;   // Final TM-score (TM_0 from USalign)
    double rmsd;     // RMSD of the final aligned region
    int    L_ali;    // Final aligned length (n_ali8 from USalign)
    double Liden;    // Number of identical residues in alignment
    std::string seqM;   // Alignment annotation string from TM-align
    std::string seqxA;  // Aligned sequence for structure 1
    std::string seqyA;  // Aligned sequence for structure 2
};

// Align backbone CA atoms of selected chains using the TM-align core
// implementation from USalign. This function is implemented in
// TMalign_wrapper.cpp and relies on the USalign headers vendored under
// src/pdb_cpp/usalign.
TMalignResult tmalign_CA(
    const Coor &coor_1,
    const Coor &coor_2,
    const std::vector<std::string> &chain_1 = {"A"},
    const std::vector<std::string> &chain_2 = {"A"},
    int mm = 0);

// Secondary structure assignment helpers (implemented in TMalign_wrapper.cpp)
std::vector<std::string> compute_SS(const Model &model, bool gap_in_seq = false);
std::vector<std::vector<std::string>> compute_SS(const Coor &coor, bool gap_in_seq = false);

#endif // PDBCPP_TMALIGN_IFACE_H
