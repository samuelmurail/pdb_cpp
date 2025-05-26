// Taken and adapted from:
// https://github.com/pylelab/USalign

// You should cite the original TMalign paper if you use this functions in your work:
// Zhang, Y., & Skolnick, J. NAR (2005). TM-align: a protein structure alignment
// algorithm based on the TM-score.
// https://doi.org/10.1093/nar/gki524

#ifndef TMALIGN_H
#define TMALIGN_H

#include <string>
#include <vector>

#include "Model.h"
#include "geom.h"
#include "Coor.h"

char sec_str(double dis13, double dis14, double dis15,
            double dis24, double dis25, double dis35)
{
    // Assign secondary structure based on distances   
    double delta=2.1;
    if (fabs(dis15-6.37)<delta && fabs(dis14-5.18)<delta && 
        fabs(dis25-5.18)<delta && fabs(dis13-5.45)<delta &&
        fabs(dis24-5.45)<delta && fabs(dis35-5.45)<delta) {
        return 'H'; // helix
    }

    delta=1.42;
    if (fabs(dis15-13  )<delta && fabs(dis14-10.4)<delta &&
        fabs(dis25-10.4)<delta && fabs(dis13-6.1 )<delta &&
        fabs(dis24-6.1 )<delta && fabs(dis35-6.1 )<delta) {
        return 'E'; //strand
    }

    if (dis15 < 8) return 'T'; //turn
    return 'C';
}


/* secondary structure assignment for protein:
 * 1->coil, 2->helix, 3->turn, 4->strand */
vector<string> compute_SS(const Model &model, bool gap_in_seq=false) {
    size_t j1, j2, j3, j4, j5, gap_num;
    double d13, d14, d15, d24, d25, d35;

    vector<int> CA_indexes = model.get_index_select("name CA");    
    size_t res_num = CA_indexes.size();

    vector<array<char, 2>> uniq_chains= model.get_uniq_chain();
    vector<array<char, 2>> chain_array = model.get_chain();
    vector<int> resid_array = model.get_resid();

    array<char, 2> old_chain = chain_array[CA_indexes[2]];
    // Get the index of the old chain in the unique chains
    auto it = find(uniq_chains.begin(), uniq_chains.end(), old_chain);
    if (it == uniq_chains.end()) {
        throw runtime_error("Chain not found in unique chains");
    }
    int chain_index = distance(uniq_chains.begin(), it);

    vector<string> seq_vec;
    int old_resid = resid_array[CA_indexes[2]];
    seq_vec.emplace_back("CC"); // Start with 'CC' to handle the first two residues

    for(size_t i=2; i < res_num - 2; i++) {
        if (chain_array[CA_indexes[i]] != old_chain) {
            // Add the two last residues of the previous chain
            if (!seq_vec[chain_index].empty())
                seq_vec[chain_index].pop_back();
            if (!seq_vec[chain_index].empty())
                seq_vec[chain_index].pop_back();
            seq_vec[chain_index] += "CC"; // Add 'CC' for the last two residues of the previous chain
            // New chain or new residue
            old_chain = chain_array[CA_indexes[i]];
            old_resid = resid_array[CA_indexes[i]];
            // Get the index of the old chain in the unique chains
            it = find(uniq_chains.begin(), uniq_chains.end(), old_chain);
            if (it == uniq_chains.end()) {
                throw runtime_error("Chain not found in unique chains");
            }
            chain_index = distance(uniq_chains.begin(), it);
            seq_vec.emplace_back("CC");
            i += 2; // Skip the next two residues since we are adding 'CC'
        }

        j1=CA_indexes[i-2];
        j2=CA_indexes[i-1];
        j3=CA_indexes[i];
        j4=CA_indexes[i+1];
        j5=CA_indexes[i+2];        
        d13=model.distance(j1, j3);
        d14=model.distance(j1, j4);
        d15=model.distance(j1, j5);
        d24=model.distance(j2, j4);
        d25=model.distance(j2, j5);
        d35=model.distance(j3, j5);
        // cout << "i: " << i << " d13: " << d13 << " d14: " << d14
        //      << " d15: " << d15 << " d24: " << d24 << " d25: " << d25
        //      << " d35: " << d35 << endl;

        if (resid_array[CA_indexes[i]] != old_resid) {
                // New residue
                if (gap_in_seq) {
                    gap_num = resid_array[CA_indexes[i]] - old_resid;
                    cout << "gap_num: " << gap_num << endl;
                    for (size_t j = 0; j < gap_num; ++j) {
                        seq_vec[chain_index] += "-"; 
                    }
                }
                old_resid = resid_array[CA_indexes[i]];
            }

        seq_vec[chain_index] += sec_str(d13, d14, d15, d24, d25, d35);
        old_resid += 1;
    } 
    seq_vec[chain_index] += "CC"; // Add 'CC' for the last two residues of the previous chain

    return seq_vec;
}

vector<vector<string>> compute_SS(const Coor &coor, bool gap_in_seq=false) {

    vector<vector<string>> ss_vec;
    
    for (size_t i = 0; i < coor.model_size(); ++i) {
        ss_vec.push_back(compute_SS(coor.get_Models(i), gap_in_seq));
    }
    return ss_vec;
}


#endif // TMALIGN_H