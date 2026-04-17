#ifndef HBOND_H
#define HBOND_H

#pragma once
#include <cmath>
#include <string>
#include <vector>
#include <array>
#include <unordered_map>
#include <limits>
#include <cstring>

#include "Model.h"
#include "geom.h"

// ─── Result type ──────────────────────────────────────────────────────────────

struct HBond {
    // Donor side
    int   donor_resid;           // uniq_resid of donor residue
    std::string donor_resname;
    std::string donor_chain;
    std::string donor_heavy_name; // e.g. "N", "OG"
    std::string donor_h_name;     // e.g. "H", "HG" (may be "virtual" if reconstructed)
    std::array<float, 3> donor_heavy_xyz;
    std::array<float, 3> donor_h_xyz;   // actual or reconstructed

    // Acceptor side
    int   acceptor_resid;
    std::string acceptor_resname;
    std::string acceptor_chain;
    std::string acceptor_name;
    std::array<float, 3> acceptor_xyz;

    // Geometry
    float dist_DA;   // D···A distance (Å)
    float dist_HA;   // H···A distance (Å)
    float angle_DHA; // D−H···A angle (degrees)
};

// ─── Protonation / acceptor tables ───────────────────────────────────────────

// Each entry: {heavy_donor_atom, hydrogen_name}
// The hydrogen_name may be empty ("") if the H must be inferred (backbone).
// Backbone N-H is handled separately at runtime using sequential residue geometry.
using DonorEntry  = std::pair<std::string, std::string>; // {heavy, H_name}
using AcceptorEntry = std::string;                        // acceptor heavy atom name

inline std::unordered_map<std::string, std::vector<DonorEntry>> build_donor_table() {
    // Amino acids — backbone N handled dynamically; only sidechain donors listed here.
    // H names follow standard PDB/CIF naming.
    std::unordered_map<std::string, std::vector<DonorEntry>> t;

    // Glycine: only backbone N
    t["GLY"] = {};
    // Alanine, Val, Leu, Ile, Pro (no sidechain donors; PRO has no backbone NH)
    t["ALA"] = {};
    t["VAL"] = {};
    t["LEU"] = {};
    t["ILE"] = {};
    t["PRO"] = {};   // N is tertiary, no H
    t["MET"] = {};
    t["PHE"] = {};
    t["TRP"] = {{"NE1", "HE1"}};
    t["TYR"] = {{"OH",  "HH"}};
    t["SER"] = {{"OG",  "HG"}};
    t["THR"] = {{"OG1", "HG1"}};
    t["CYS"] = {{"SG",  "HG"}};
    t["ASN"] = {{"ND2", "HD21"}, {"ND2", "HD22"}};
    t["GLN"] = {{"NE2", "HE21"}, {"NE2", "HE22"}};
    t["HIS"] = {{"ND1", "HD1"},  {"NE2", "HE2"}};  // either or both depending on protonation
    t["HSD"] = {{"ND1", "HD1"}};
    t["HSE"] = {{"NE2", "HE2"}};
    t["HSP"] = {{"ND1", "HD1"}, {"NE2", "HE2"}};
    t["HID"] = {{"ND1", "HD1"}};
    t["HIE"] = {{"NE2", "HE2"}};
    t["HIP"] = {{"ND1", "HD1"}, {"NE2", "HE2"}};
    t["LYS"] = {{"NZ", "HZ1"}, {"NZ", "HZ2"}, {"NZ", "HZ3"}};
    t["ARG"] = {{"NH1", "HH11"}, {"NH1", "HH12"}, {"NH2", "HH21"}, {"NH2", "HH22"}, {"NE", "HE"}};
    t["ASP"] = {};   // deprotonated; protonated form:
    t["ASH"] = {{"OD1", "HD1"}, {"OD2", "HD2"}};
    t["ASPP"]= {{"OD1", "HD1"}, {"OD2", "HD2"}};
    t["GLU"] = {};
    t["GLH"] = {{"OE1", "HE1"}, {"OE2", "HE2"}};
    t["GLUP"]= {{"OE1", "HE1"}, {"OE2", "HE2"}};
    t["SEC"] = {{"SE",  "HSE"}};
    // D-amino acids (same donors, different residue names)
    t["DAL"] = {};
    t["DAR"] = t["ARG"];
    t["DSG"] = t["ASN"];
    t["DAS"] = {};
    t["DCY"] = t["CYS"];
    t["DGN"] = t["GLN"];
    t["DGL"] = {};
    t["DHI"] = t["HIS"];
    t["DIL"] = {};
    t["DLE"] = {};
    t["DLY"] = t["LYS"];
    t["DME"] = {};
    t["MED"] = {};
    t["DPH"] = {};
    t["DPN"] = {};
    t["DPR"] = {};
    t["DSE"] = t["SER"];
    t["DSN"] = t["SER"];
    t["DTH"] = t["THR"];
    t["DTR"] = t["TRP"];
    t["DTY"] = t["TYR"];
    t["DVA"] = {};
    // RNA
    t["A"]  = {{"N6", "H61"}, {"N6", "H62"}};
    t["U"]  = {{"N3", "H3"}};
    t["G"]  = {{"N1", "H1"}, {"N2", "H21"}, {"N2", "H22"}};
    t["C"]  = {{"N4", "H41"}, {"N4", "H42"}};
    // DNA
    t["DA"] = t["A"];
    t["DT"] = {{"N3", "H3"}};
    t["DC"] = t["C"];
    t["DG"] = t["G"];

    return t;
}

inline std::unordered_map<std::string, std::vector<AcceptorEntry>> build_acceptor_table() {
    std::unordered_map<std::string, std::vector<AcceptorEntry>> t;

    // All amino acids: backbone O is always an acceptor
    const std::vector<std::string> all_aa = {
        "GLY","ALA","VAL","LEU","ILE","PRO","MET","PHE","TRP","TYR",
        "SER","THR","CYS","ASN","GLN","HIS","HSD","HSE","HSP","HID","HIE","HIP",
        "LYS","ARG","ASP","ASH","ASPP","GLU","GLH","GLUP","SEC",
        "DAL","DAR","DSG","DAS","DCY","DGN","DGL","DHI","DIL","DLE",
        "DLY","DME","MED","DPH","DPN","DPR","DSE","DSN","DTH","DTR","DTY","DVA"
    };
    for (const auto &aa : all_aa) {
        t[aa] = {"O"};  // backbone carbonyl O
    }

    // Sidechain acceptors
    t["SER"].push_back("OG");
    t["DSE"].push_back("OG");
    t["DSN"].push_back("OG");
    t["THR"].push_back("OG1");
    t["DTH"].push_back("OG1");
    t["TYR"].push_back("OH");
    t["DTY"].push_back("OH");
    t["CYS"].push_back("SG");
    t["DCY"].push_back("SG");
    t["MET"].push_back("SD");
    t["DME"].push_back("SD");
    t["MED"].push_back("SD");
    t["ASP"].push_back("OD1"); t["ASP"].push_back("OD2");
    t["ASH"].push_back("OD1"); t["ASH"].push_back("OD2");
    t["ASPP"].push_back("OD1"); t["ASPP"].push_back("OD2");
    t["DAS"].push_back("OD1"); t["DAS"].push_back("OD2");
    t["GLU"].push_back("OE1"); t["GLU"].push_back("OE2");
    t["GLH"].push_back("OE1"); t["GLH"].push_back("OE2");
    t["GLUP"].push_back("OE1"); t["GLUP"].push_back("OE2");
    t["DGL"].push_back("OE1"); t["DGL"].push_back("OE2");
    t["ASN"].push_back("OD1");
    t["DSG"].push_back("OD1");
    t["GLN"].push_back("OE1");
    t["DGN"].push_back("OE1");
    t["HIS"].push_back("ND1"); t["HIS"].push_back("NE2");
    t["HSD"].push_back("NE2");
    t["HSE"].push_back("ND1");
    t["HID"].push_back("NE2");
    t["HIE"].push_back("ND1");
    t["DHI"].push_back("ND1"); t["DHI"].push_back("NE2");
    t["SEC"].push_back("SE");
    // RNA / DNA base acceptors
    t["A"]  = {"N1", "N3", "N7"};
    t["U"]  = {"O2", "O4"};
    t["G"]  = {"N3", "N7", "O6"};
    t["C"]  = {"N3", "O2"};
    t["DA"] = t["A"];
    t["DT"] = {"O2", "O4"};
    t["DC"] = t["C"];
    t["DG"] = t["G"];

    // Backbone / sugar acceptors.
    // Phosphate oxygens are a major class of protein–nucleic H-bond acceptors.
    const std::vector<std::string> nucleic_all = {"A", "U", "G", "C", "DA", "DT", "DC", "DG"};
    const std::vector<std::string> backbone_acceptors = {"OP1", "OP2", "O1P", "O2P", "O3'", "O4'", "O5'"};
    for (const auto &nt : nucleic_all) {
        auto &acc = t[nt];
        acc.insert(acc.end(), backbone_acceptors.begin(), backbone_acceptors.end());
    }

    // RNA 2'-OH oxygen can accept as well.
    t["A"].push_back("O2'");
    t["U"].push_back("O2'");
    t["G"].push_back("O2'");
    t["C"].push_back("O2'");

    return t;
}

// ─── Helpers ─────────────────────────────────────────────────────────────────

inline std::string char5_to_str(const std::array<char, 5> &a) {
    return std::string(a.data(), strnlen(a.data(), 5));
}
inline std::string char2_to_str(const std::array<char, 2> &a) {
    if (a[0] == ' ' || a[0] == '\0') return "";
    if (a[1] == '\0') return std::string(1, a[0]);
    return std::string(a.data(), 2);
}

inline float dot3(const float *a, const float *b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
inline float norm3(const float *a) {
    return std::sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
}
inline float angle_DHA(const std::array<float,3> &D, const std::array<float,3> &H, const std::array<float,3> &A) {
    // Angle at H: D-H···A
    float DH[3] = {D[0]-H[0], D[1]-H[1], D[2]-H[2]};
    float AH[3] = {A[0]-H[0], A[1]-H[1], A[2]-H[2]};
    float n_DH = norm3(DH);
    float n_AH = norm3(AH);
    if (n_DH < 1e-6f || n_AH < 1e-6f) return 0.0f;
    float cosA = dot3(DH, AH) / (n_DH * n_AH);
    if (cosA >  1.0f) cosA =  1.0f;
    if (cosA < -1.0f) cosA = -1.0f;
    return std::acos(cosA) * (180.0f / 3.14159265358979323846f);
}

// ─── Per-residue index helpers ───────────────────────────────────────────────

// Returns {x, y, z} for the first atom in `model` at uniq_resid==uid with name==aname,
// or {NaN,NaN,NaN} if not found.
inline std::array<float,3> find_atom_xyz(
        const Model &model, int uid, const std::string &aname) {
    constexpr float NaN = std::numeric_limits<float>::quiet_NaN();
    const auto &uids  = model.get_uniqresid();
    const auto &names = model.get_name();
    const auto &xs    = model.get_x();
    const auto &ys    = model.get_y();
    const auto &zs    = model.get_z();
    for (size_t i = 0; i < uids.size(); ++i) {
        if (uids[i] == uid && char5_to_str(names[i]) == aname) {
            return {xs[i], ys[i], zs[i]};
        }
    }
    return {NaN, NaN, NaN};
}

// Reconstruct backbone N-H virtual hydrogen for residue at `uid`,
// using C of residue at `prev_uid`. Returns NaN array if precursors missing.
template<typename FindFn>
inline std::array<float,3> reconstruct_backbone_NH(
        FindFn &&find, int uid, int prev_uid) {
    constexpr float NaN = std::numeric_limits<float>::quiet_NaN();
    auto N  = find(uid, "N");
    auto CA = find(uid, "CA");
    auto C_prev = find(prev_uid, "C");
    if (std::isnan(N[0]) || std::isnan(CA[0]) || std::isnan(C_prev[0]))
        return {NaN, NaN, NaN};

    // Unit vector N → C(i-1)
    float CN[3] = {N[0]-C_prev[0], N[1]-C_prev[1], N[2]-C_prev[2]};
    float n = norm3(CN);
    if (n < 1e-6f) return {NaN, NaN, NaN};
    CN[0] /= n; CN[1] /= n; CN[2] /= n;

    // Unit vector N → CA
    float CaN[3] = {N[0]-CA[0], N[1]-CA[1], N[2]-CA[2]};
    float m = norm3(CaN);
    if (m < 1e-6f) return {NaN, NaN, NaN};
    CaN[0] /= m; CaN[1] /= m; CaN[2] /= m;

    // Bisector → H direction, 1.01 Å
    float H[3] = {CN[0]+CaN[0], CN[1]+CaN[1], CN[2]+CaN[2]};
    float hn = norm3(H);
    if (hn < 1e-6f) return {NaN, NaN, NaN};
    return {N[0] + 1.01f*H[0]/hn,
            N[1] + 1.01f*H[1]/hn,
            N[2] + 1.01f*H[2]/hn};
}

// Reconstruct a sidechain H from the known heavy atom and one bonded neighbour,
// using sp3-like tetrahedral geometry (bond angle ~109.5°, bond length 1.0 Å).
// This is a minimal single-neighbour version for OH/NH2/SH/NH groups.
// For atoms with a single known bonded neighbour B–D–H, the H is placed along
// the extension of B→D at 1.0 Å.
inline std::array<float,3> reconstruct_sidechain_H(
        const std::array<float,3> &D, const std::array<float,3> &B,
        float bond_len = 1.0f) {
    float v[3] = {D[0]-B[0], D[1]-B[1], D[2]-B[2]};
    float n = norm3(v);
    if (n < 1e-6f) {
        constexpr float NaN = std::numeric_limits<float>::quiet_NaN();
        return {NaN, NaN, NaN};
    }
    return {D[0] + bond_len*v[0]/n,
            D[1] + bond_len*v[1]/n,
            D[2] + bond_len*v[2]/n};
}

// ─── Donor list builder ───────────────────────────────────────────────────────

// Returns list of (heavy_xyz, H_xyz, heavy_name, H_name) for all donors in
// residue uid in model, given donor table.
// Backbone N-H is added if prev_uid >= 0 and residue is not PRO.
struct DonorAtom {
    std::array<float,3> heavy_xyz;
    std::array<float,3> h_xyz;
    std::string heavy_name;
    std::string h_name;
};

template<typename FindFn>
inline std::vector<DonorAtom> get_donors_for_residue(
        FindFn &&find,
        int uid,
        const std::string &resname,
        int prev_uid,
        const std::unordered_map<std::string, std::vector<DonorEntry>> &donor_table) {

    constexpr float NaN = std::numeric_limits<float>::quiet_NaN();
    std::vector<DonorAtom> donors;

    // 1. Backbone N-H (all amino acids except PRO; not first residue)
    const bool is_aa = !(resname == "DA" || resname == "DT" || resname == "DC" || resname == "DG" ||
                         resname == "A"  || resname == "U"  || resname == "G"  || resname == "C");
    if (is_aa && resname != "PRO" && resname != "DPR" && prev_uid >= 0) {
        auto N_xyz = find(uid, "N");
        if (!std::isnan(N_xyz[0])) {
            // Try to find H directly in coordinates first
            auto H_xyz = find(uid, "H");
            if (std::isnan(H_xyz[0]))
                H_xyz = reconstruct_backbone_NH(find, uid, prev_uid);
            if (!std::isnan(H_xyz[0]))
                donors.push_back({N_xyz, H_xyz, "N", "H"});
        }
    }

    // 2. Sidechain donors from table
    auto it = donor_table.find(resname);
    if (it == donor_table.end()) return donors;
    for (const auto &entry : it->second) {
        const std::string &heavy = entry.first;
        const std::string &hname = entry.second;
        auto D_xyz = find(uid, heavy);
        if (std::isnan(D_xyz[0])) continue;

        // Try to read H from model coords
        std::array<float,3> H_xyz = {NaN, NaN, NaN};
        if (!hname.empty())
            H_xyz = find(uid, hname);

        if (std::isnan(H_xyz[0])) {
            // Geometrically reconstruct: need one bonded neighbour
            // Standard bonded neighbours by heavy atom role:
            std::string neighbour;
            if (heavy == "N" || heavy == "ND1" || heavy == "NE2" || heavy == "NE" || heavy == "NZ" ||
                heavy == "NH1" || heavy == "NH2" || heavy == "ND2" || heavy == "NE1")
                neighbour = "CA";  // fallback: connect N to CA
            else if (heavy == "OG" || heavy == "OG1")
                neighbour = "CB";
            else if (heavy == "OH")
                neighbour = "CZ";
            else if (heavy == "SG")
                neighbour = "CB";
            else if (heavy == "OD1" || heavy == "OD2")
                neighbour = "CG";
            else if (heavy == "OE1" || heavy == "OE2")
                neighbour = "CD";
            else if (heavy == "N1" || heavy == "N2" || heavy == "N3" || heavy == "N4" || heavy == "N6")
                neighbour = "C1'"; // nucleic acid
            else
                neighbour = "";

            if (!neighbour.empty()) {
                auto B_xyz = find(uid, neighbour);
                if (!std::isnan(B_xyz[0]))
                    H_xyz = reconstruct_sidechain_H(D_xyz, B_xyz);
            }
        }
        if (!std::isnan(H_xyz[0]))
            donors.push_back({D_xyz, H_xyz, heavy, hname.empty() ? "H" : hname});
    }
    return donors;
}

// ─── Main function ────────────────────────────────────────────────────────────

/// Compute H-bonds between atoms in `donor_model` (donors) and `acceptor_model` (acceptors).
/// Both models are query Model objects (typically from Coor::select_atoms).
/// `full_model` is the full frame, used to find prev-residue C for backbone H reconstruction.
/// Criteria (Baker & Hubbard):
///   D···A distance  < dist_DA_cutoff  (default 3.5 Å)
///   H···A distance  < dist_HA_cutoff  (default 2.5 Å)
///   D−H···A angle   > angle_cutoff    (default 90°)
inline std::vector<HBond> compute_hbonds(
        const Model &donor_model,
        const Model &acceptor_model,
        const Model &full_model,
        float dist_DA_cutoff = 3.5f,
        float dist_HA_cutoff = 2.5f,
        float angle_cutoff   = 90.0f) {

    static const auto DONOR_TABLE    = build_donor_table();
    static const auto ACCEPTOR_TABLE = build_acceptor_table();

    constexpr float NaN = std::numeric_limits<float>::quiet_NaN();
    std::vector<HBond> result;

    // ── collect unique residues from donor_model ──────────────────────────
    const auto &d_uids    = donor_model.get_uniqresid();
    const auto &d_resnames= donor_model.get_resname();
    const auto &d_chains  = donor_model.get_chain();

    // build ordered unique-resid list from donor_model
    std::vector<int> donor_uids_ordered;
    {
        std::unordered_map<int, bool> seen;
        for (int uid : d_uids)
            if (!seen[uid]) { seen[uid]=true; donor_uids_ordered.push_back(uid); }
    }

    // map uid → resname / chain (take first occurrence)
    std::unordered_map<int, std::string> uid_to_resname;
    std::unordered_map<int, std::string> uid_to_chain;
    for (size_t i = 0; i < d_uids.size(); ++i) {
        int uid = d_uids[i];
        if (uid_to_resname.find(uid) == uid_to_resname.end()) {
            uid_to_resname[uid] = char5_to_str(d_resnames[i]);
            uid_to_chain[uid]   = char2_to_str(d_chains[i]);
        }
    }

    // build prev_uid map from full_model (ordered unique residue list)
    const auto &fm_uids = full_model.get_uniqresid();
    std::vector<int> full_uids_ordered;
    {
        std::unordered_map<int, bool> seen;
        for (int uid : fm_uids)
            if (!seen[uid]) { seen[uid]=true; full_uids_ordered.push_back(uid); }
    }
    std::unordered_map<int, int> uid_to_prev;
    for (size_t i = 0; i < full_uids_ordered.size(); ++i)
        uid_to_prev[full_uids_ordered[i]] = (i == 0) ? -1 : full_uids_ordered[i-1];

    // ── pre-build (uid, atom_name) → xyz lookup from full_model ──────────
    // Replaces the O(N_atoms) linear scan inside find_atom_xyz with O(1) lookup.
    struct UidName {
        int uid; std::string name;
        bool operator==(const UidName &o) const { return uid == o.uid && name == o.name; }
    };
    struct UidNameHash {
        std::size_t operator()(const UidName &k) const noexcept {
            std::size_t h = std::hash<int>{}(k.uid);
            for (char c : k.name) h = h * 31 + (unsigned char)c;
            return h;
        }
    };
    std::unordered_map<UidName, std::array<float,3>, UidNameHash> atom_xyz_map;
    {
        const auto &fm_names = full_model.get_name();
        const auto &fm_xs    = full_model.get_x();
        const auto &fm_ys    = full_model.get_y();
        const auto &fm_zs    = full_model.get_z();
        atom_xyz_map.reserve(fm_uids.size());
        for (size_t i = 0; i < fm_uids.size(); ++i) {
            UidName key{fm_uids[i], char5_to_str(fm_names[i])};
            atom_xyz_map.emplace(key, std::array<float,3>{fm_xs[i], fm_ys[i], fm_zs[i]});
        }
    }

    // Lambda replacing find_atom_xyz — O(1) instead of O(N_atoms)
    auto fast_find = [&](int uid, const std::string &aname) -> std::array<float,3> {
        constexpr float NaN = std::numeric_limits<float>::quiet_NaN();
        auto it = atom_xyz_map.find({uid, aname});
        return (it != atom_xyz_map.end()) ? it->second
                                          : std::array<float,3>{NaN, NaN, NaN};
    };

    // ── collect acceptors from acceptor_model ─────────────────────────────
    const auto &a_uids    = acceptor_model.get_uniqresid();
    const auto &a_resnames= acceptor_model.get_resname();
    const auto &a_names   = acceptor_model.get_name();
    const auto &a_chains  = acceptor_model.get_chain();
    const auto &a_xs      = acceptor_model.get_x();
    const auto &a_ys      = acceptor_model.get_y();
    const auto &a_zs      = acceptor_model.get_z();

    struct AcceptorAtom {
        int uid;
        std::string resname;
        std::string chain;
        std::string name;
        std::array<float,3> xyz;
    };
    std::vector<AcceptorAtom> acceptors;
    for (size_t i = 0; i < a_uids.size(); ++i) {
        std::string rn = char5_to_str(a_resnames[i]);
        std::string an = char5_to_str(a_names[i]);
        auto it = ACCEPTOR_TABLE.find(rn);
        if (it == ACCEPTOR_TABLE.end()) continue;
        for (const auto &acc : it->second) {
            if (an == acc) {
                acceptors.push_back({a_uids[i], rn, char2_to_str(a_chains[i]),
                                     an, {a_xs[i], a_ys[i], a_zs[i]}});
            }
        }
    }

    // ── build cell list for acceptors (voxel size = dist_DA_cutoff) ───────
    // This turns the O(N_donors × N_acceptors) search into O(N · k) where
    // k is the average number of acceptors in the 3×3×3 neighbourhood.
    const float cell = dist_DA_cutoff;   // one voxel side length
    using CellKey = std::tuple<int,int,int>;
    struct CellKeyHash {
        std::size_t operator()(const CellKey &k) const noexcept {
            // FNV-like combine
            std::size_t h = (std::size_t)(std::get<0>(k) * 73856093)
                          ^ (std::size_t)(std::get<1>(k) * 19349663)
                          ^ (std::size_t)(std::get<2>(k) * 83492791);
            return h;
        }
    };
    std::unordered_map<CellKey, std::vector<size_t>, CellKeyHash> cell_map;
    cell_map.reserve(acceptors.size() * 2);
    for (size_t i = 0; i < acceptors.size(); ++i) {
        const auto &xyz = acceptors[i].xyz;
        int cx = (int)std::floor(xyz[0] / cell);
        int cy = (int)std::floor(xyz[1] / cell);
        int cz = (int)std::floor(xyz[2] / cell);
        cell_map[{cx, cy, cz}].push_back(i);
    }

    // ── for each donor residue, enumerate donor atoms, check against acceptors ──
    for (int d_uid : donor_uids_ordered) {
        auto rn_it = uid_to_resname.find(d_uid);
        if (rn_it == uid_to_resname.end()) continue;
        const std::string &resname = rn_it->second;
        const std::string &chain   = uid_to_chain[d_uid];
        int prev_uid = uid_to_prev.count(d_uid) ? uid_to_prev[d_uid] : -1;

        auto donor_atoms = get_donors_for_residue(
            fast_find, d_uid, resname, prev_uid, DONOR_TABLE);

        for (const auto &don : donor_atoms) {
            const auto &D = don.heavy_xyz;
            const auto &H = don.h_xyz;

            // Query the 3×3×3 neighbourhood of donor D in the cell list
            int dcx = (int)std::floor(D[0] / cell);
            int dcy = (int)std::floor(D[1] / cell);
            int dcz = (int)std::floor(D[2] / cell);

            for (int ox = -1; ox <= 1; ++ox)
            for (int oy = -1; oy <= 1; ++oy)
            for (int oz = -1; oz <= 1; ++oz) {
                auto cit = cell_map.find({dcx+ox, dcy+oy, dcz+oz});
                if (cit == cell_map.end()) continue;
                for (size_t ai : cit->second) {
                    const auto &acc = acceptors[ai];
                // Skip self-residue
                if (acc.uid == d_uid) continue;

                const auto &A = acc.xyz;

                // Quick D···A pre-filter
                float dx = D[0]-A[0], dy = D[1]-A[1], dz = D[2]-A[2];
                float dDA = std::sqrt(dx*dx + dy*dy + dz*dz);
                if (dDA > dist_DA_cutoff) continue;

                // H···A distance
                float hx = H[0]-A[0], hy = H[1]-A[1], hz = H[2]-A[2];
                float dHA = std::sqrt(hx*hx + hy*hy + hz*hz);
                if (dHA > dist_HA_cutoff) continue;

                // D−H···A angle
                float ang = angle_DHA(D, H, A);
                if (ang < angle_cutoff) continue;

                HBond hb;
                hb.donor_resid      = d_uid;
                hb.donor_resname    = resname;
                hb.donor_chain      = chain;
                hb.donor_heavy_name = don.heavy_name;
                hb.donor_h_name     = don.h_name;
                hb.donor_heavy_xyz  = D;
                hb.donor_h_xyz      = H;
                hb.acceptor_resid   = acc.uid;
                hb.acceptor_resname = acc.resname;
                hb.acceptor_chain   = acc.chain;
                hb.acceptor_name    = acc.name;
                hb.acceptor_xyz     = A;
                hb.dist_DA          = dDA;
                hb.dist_HA          = dHA;
                hb.angle_DHA        = ang;
                result.push_back(hb);
                } // for ai in cell
            } // for oz (closes ox/oy/oz block)
        } // for don
    } // for d_uid
    return result;
}

#endif // HBOND_H
