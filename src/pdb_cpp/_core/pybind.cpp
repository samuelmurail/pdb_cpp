#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include "Model.h"
#include "Coor.h"
#include "TMalign_iface.h"
#include "align.h"
#include "geom.h"
#include "seq_align.h"
#include "format/encode.h"
#include "hbond.h"

namespace py = pybind11;

static std::string resolve_matrix_file(const std::string &matrix_file) {
    if (!matrix_file.empty()) {
        return matrix_file;
    }
    py::module_ resources = py::module_::import("importlib.resources");
    py::object matrix_path = resources.attr("files")("pdb_cpp.data").attr("joinpath")("blosum62.txt");
    return py::str(matrix_path);
}

PYBIND11_MODULE(core, m) {
    py::class_<Model>(m, "Model")
        .def(py::init<>())
        .def("clear", &Model::clear)
        .def("size", &Model::size)
        .def("addAtom", &Model::addAtom)
        .def("select_atoms", &Model::select_atoms)
        

        .def("get_x", &Model::get_x)
        .def("get_y", &Model::get_y)
        .def("get_z", &Model::get_z)
        .def("get_name", &Model::get_name)
        .def("get_resname", &Model::get_resname)
        .def("get_resid", &Model::get_resid)
        .def("get_chain", &Model::get_chain)
        .def("get_occ", &Model::get_occ)
        .def("get_beta", &Model::get_beta)
        .def("get_alterloc", &Model::get_alterloc)
        .def("get_insertres", &Model::get_insertres)
        .def("get_elem", &Model::get_elem)
        .def("get_num", &Model::get_num)
        .def("get_field", &Model::get_field)
        .def("get_uniqresid", &Model::get_uniqresid)
        .def("get_centroid", static_cast<std::array<float, 3> (Model::*)() const>(&Model::get_centroid),
             "Calculate centroid of all atoms in the model")
        .def("get_centroid", static_cast<std::array<float, 3> (Model::*)(const std::vector<int>&) const>(&Model::get_centroid),
             py::arg("indices"), "Calculate centroid of atoms at specified indices")
        // Add more methods as needed
        ;
    py::class_<Coor>(m, "Coor")
        .def_property("active_model", &Coor::get_active_model, &Coor::set_active_model)
        .def(py::init<>())
        .def(py::init(
            [](const std::string &filename, const std::string &format) {
                Coor c;
                c.read(filename, format);
                return c;
            }),
            py::arg("filename"), py::arg("format") = "",
            "Construct a Coor from a file; format can be 'pdb', 'cif', 'pqr', or 'gro' (default: infer from extension)")
        .def("read", &Coor::read,
            py::arg("filename"), py::arg("format") = "",
            "Read a structure file; format can be 'pdb', 'cif', 'pqr', or 'gro' (default: infer from extension)")
        .def("write", &Coor::write)
        .def("clear", &Coor::clear)
        .def("size", &Coor::size)
        .def("add_Model", &Coor::add_Model)
        .def("get_Models", &Coor::get_Models)
        .def("get_all_Models", &Coor::get_all_Models)
        .def("model_size", &Coor::model_size)
        .def_readwrite("conect", &Coor::conect)
        .def("select_atoms", &Coor::select_atoms, 
            py::arg("selection"), py::arg("frame") = 0, // Specify default value for `frame`
            "Select atoms based on a selection string and an optional frame index")
        .def("get_index_select", &Coor::get_index_select,
            py::arg("selection"), py::arg("frame") = 0, // Specify default value for `frame`
            "Get indices of atoms based on a selection string and an optional frame index")
        .def("select_bool_index", &Coor::select_bool_index)
        .def("get_uniq_chain", &Coor::get_uniq_chain)
        .def("get_uniq_chain_str", &Coor::get_uniq_chain_str)
        .def("get_aa_sequences", &Coor::get_aa_sequences, 
            py::arg("gap_in_seq") = true, py::arg("frame") = 0, // Specify default values for `gap_in_seq` and `frame`
            "Get the amino acid sequence, optionally including gaps and specifying a frame index")
        .def("get_aa_sequences_dl", &Coor::get_aa_sequences_dl,
            py::arg("gap_in_seq") = true, py::arg("frame") = 0,
            "Get the amino acid sequence with D-residues encoded as lowercase")
        .def("get_aa_seq", &Coor::get_aa_seq,
            py::arg("gap_in_seq") = true, py::arg("frame") = 0,
            "Get amino-acid sequences per chain")
        .def("get_aa_DL_seq", &Coor::get_aa_DL_seq,
            py::arg("gap_in_seq") = true, py::arg("frame") = 0,
            "Get amino-acid sequences with D-residues encoded as lowercase")
        .def("get_aa_na_seq", &Coor::get_aa_na_seq,
            py::arg("gap_in_seq") = true, py::arg("frame") = 0,
            "Get amino-acid and nucleic-acid sequences per chain")
        .def("remove_incomplete_backbone_residues", &Coor::remove_incomplete_backbone_residues,
            py::arg("back_atom") = std::vector<std::string>{"CA", "C", "N", "O"},
            "Remove residues with incomplete backbone atoms")
        // Add more methods as needed
        ;

    m.def("distance_matrix",
        [](py::array_t<float, py::array::c_style | py::array::forcecast> a,
           py::array_t<float, py::array::c_style | py::array::forcecast> b) {
            if (a.ndim() != 2 || b.ndim() != 2 || a.shape(1) != 3 || b.shape(1) != 3) {
                throw std::runtime_error("Inputs must be 2D arrays with shape (N, 3) and (M, 3)");
            }
            const auto n = static_cast<size_t>(a.shape(0));
            const auto m = static_cast<size_t>(b.shape(0));
            py::array_t<float> out({n, m});
            const float *a_ptr = static_cast<const float *>(a.data());
            const float *b_ptr = static_cast<const float *>(b.data());
            float *out_ptr = static_cast<float *>(out.mutable_data());
            distance_matrix(a_ptr, n, b_ptr, m, out_ptr);
            return out;
        },
        py::arg("xyz_a"),
        py::arg("xyz_b"),
        "Compute a pairwise distance matrix between two coordinate sets");

    m.def("compute_dihedrals",
        [](py::array_t<float, py::array::c_style | py::array::forcecast> pts) {
            if (pts.ndim() != 2 || pts.shape(1) != 3)
                throw std::runtime_error("Input must be a 2D array with shape (N, 3)");
            const auto n = static_cast<size_t>(pts.shape(0));
            if (n < 4) return py::array_t<float>(py::array::ShapeContainer{0});
            py::array_t<float> out(py::array::ShapeContainer{n - 3});
            compute_dihedrals(pts.data(), n, out.mutable_data());
            return out;
        },
        py::arg("pts"),
        R"doc(
Compute all consecutive dihedral angles from an ordered (N, 3) float array.

Parameters
----------
pts : ndarray, shape (N, 3)
    Ordered 3-D coordinates (e.g. consecutive CA positions).

Returns
-------
ndarray, shape (N-3,)
    Dihedral angles in degrees.  Returns an empty array when N < 4.
)doc");
    // Bind the TMAlign secondary-structure function
    m.def("compute_SS",
        static_cast<std::vector<std::vector<std::string>>(*)(const Coor&, bool)>(&compute_SS),
        py::arg("coor"),
        py::arg("gap_in_seq") = false,
        "Compute secondary structure for all models in a Coor object");

    // Bind the TM-align structural alignment result
    py::class_<TMalignResult>(m, "TMalignResult")
        .def_readonly("TM1", &TMalignResult::TM1)
        .def_readonly("TM2", &TMalignResult::TM2)
        .def_readonly("TM_ali", &TMalignResult::TM_ali)
        .def_readonly("rmsd", &TMalignResult::rmsd)
        .def_readonly("L_ali", &TMalignResult::L_ali)
        .def_readonly("Liden", &TMalignResult::Liden)
        .def_readonly("seqM", &TMalignResult::seqM)
        .def_readonly("seqxA", &TMalignResult::seqxA)
        .def_readonly("seqyA", &TMalignResult::seqyA);

    // Bind the TM-align based CA alignment helper
    m.def("tmalign_ca",
        &tmalign_CA,
        py::arg("coor_1"),
        py::arg("coor_2"),
        py::arg("chain_1") = std::vector<std::string>{"A"},
        py::arg("chain_2") = std::vector<std::string>{"A"},
        py::arg("mm") = 0,
        "Align CA atoms of selected chains using the TM-align core from USalign");

    // Bind the Alignment structure
    py::class_<Alignment_cpp>(m, "Alignment_cpp")
        .def(py::init<>())  // Default constructor
        .def_readwrite("seq1", &Alignment_cpp::seq1, "Aligned sequence 1")
        .def_readwrite("seq2", &Alignment_cpp::seq2, "Aligned sequence 2")
        .def_readwrite("score", &Alignment_cpp::score, "Alignment score");

    // Bind the align function
    m.def("sequence_align",
        [](const std::string &seq1,
           const std::string &seq2,
           const std::string &matrix_file,
           int gap_cost,
           int gap_ext) {
            return sequence_align(seq1, seq2, resolve_matrix_file(matrix_file), gap_cost, gap_ext);
        },
        py::arg("seq1"),
        py::arg("seq2"),
        py::arg("matrix_file") = "",
        py::arg("GAP_COST") = -11,
        py::arg("GAP_EXT") = -1,
        py::return_value_policy::take_ownership,  // Python takes ownership of the returned pointer
        "Align two sequences using a substitution matrix and gap penalties");

    // Bind the get_common_atoms function
    m.def("get_common_atoms",
        [](const Coor &coor_1,
           const Coor &coor_2,
           const std::vector<std::string> &chain_1,
           const std::vector<std::string> &chain_2,
           const std::vector<std::string> &back_names,
           const std::string &matrix_file) {
            return get_common_atoms(
                coor_1,
                coor_2,
                chain_1,
                chain_2,
                back_names,
                resolve_matrix_file(matrix_file)
            );
        },
        py::arg("coor_1"),
        py::arg("coor_2"),
        py::arg("chain_1") = std::vector<std::string>{"A"},
        py::arg("chain_2") = std::vector<std::string>{"A"},
        py::arg("back_names") = std::vector<std::string>{"C", "N", "O", "CA"},
        py::arg("matrix_file") = "",
        "Get common atoms between two Coor objects based on sequence alignment");

    // Bind the coor_align function
    m.def("coor_align",
        &coor_align,
        py::arg("coor_1"),
        py::arg("coor_2"),
        py::arg("index_1"),
        py::arg("index_2"),
        py::arg("frame_ref") = 0,
        "Align two coordinate structures using quaternion-based rotation");

    m.def("rmsd",
        [](const Coor &coor_1,
           const Coor &coor_2,
           const std::vector<int> &index_1,
           const std::vector<int> &index_2,
           int frame_ref) {
            if (frame_ref < 0 || static_cast<size_t>(frame_ref) >= coor_2.model_size()) {
                throw std::runtime_error(
                    "Reference frame index is larger than the number of frames in the reference structure"
                );
            }
            if (index_1.empty() || index_2.empty()) {
                throw std::runtime_error("No atoms selected for RMSD calculation");
            }
            if (index_1.size() != index_2.size()) {
                throw std::runtime_error("Index lists must have the same length");
            }

            const Model ref_model = coor_2.get_Models(frame_ref);
            std::vector<float> rmsd_values;
            rmsd_values.reserve(coor_1.model_size());
            for (size_t model_index = 0; model_index < coor_1.model_size(); ++model_index) {
                const Model model = coor_1.get_Models(static_cast<int>(model_index));
                rmsd_values.push_back(::rmsd(model, ref_model, index_1, index_2));
            }
            return rmsd_values;
        },
        py::arg("coor_1"),
        py::arg("coor_2"),
        py::arg("index_1"),
        py::arg("index_2"),
        py::arg("frame_ref") = 0,
        "Compute RMSD values for all models in coor_1 against one reference model in coor_2.");
    
    m.def("align_seq_based",
        [](Coor &coor_1,
           const Coor &coor_2,
           const std::vector<std::string> &chain_1,
           const std::vector<std::string> &chain_2,
           const std::vector<std::string> &back_names,
           const std::string &matrix_file,
           int frame_ref) {
            return align_seq_based(
                coor_1,
                coor_2,
                chain_1,
                chain_2,
                back_names,
                resolve_matrix_file(matrix_file),
                frame_ref
            );
        },
        py::arg("coor_1"),
        py::arg("coor_2"),
        py::arg("chain_1") = std::vector<std::string>{"A"},
        py::arg("chain_2") = std::vector<std::string>{"A"},
        py::arg("back_names") = std::vector<std::string>{"C", "N", "O", "CA"},
        py::arg("matrix_file") = "",
        py::arg("frame_ref") = 0,
        "Align two coordinate structures using sequence based alignment");

    m.def("align_chain_permutation",
        [](const Coor &coor_1,
           const Coor &coor_2,
           const std::vector<std::string> &back_names,
           const std::string &matrix_file,
           int frame_ref) {
            return align_chain_permutation(
                coor_1,
                coor_2,
                back_names,
                resolve_matrix_file(matrix_file),
                frame_ref
            );
        },
        py::arg("coor_1"),
        py::arg("coor_2"),
        py::arg("back_names") = std::vector<std::string>{"C", "N", "O", "CA"},
        py::arg("matrix_file") = "",
        py::arg("frame_ref") = 0,
        "Align structures by permuting chain order and selecting the best RMSD.");

    m.def("hy36encode",
        &hy36encode,
        py::arg("width"),
        py::arg("value"),
        "Encode a number using hybrid-36 with fixed width.");
    m.def("hy36decode",
        &hy36decode,
        py::arg("width"),
        py::arg("value"),
        "Decode a hybrid-36 string with fixed width.");

    // HBond result struct
    py::class_<HBond>(m, "HBond")
        .def_readonly("donor_resid",       &HBond::donor_resid)
        .def_readonly("donor_resname",     &HBond::donor_resname)
        .def_readonly("donor_chain",       &HBond::donor_chain)
        .def_readonly("donor_heavy_name",  &HBond::donor_heavy_name)
        .def_readonly("donor_h_name",      &HBond::donor_h_name)
        .def_readonly("donor_heavy_xyz",   &HBond::donor_heavy_xyz)
        .def_readonly("donor_h_xyz",       &HBond::donor_h_xyz)
        .def_readonly("acceptor_resid",    &HBond::acceptor_resid)
        .def_readonly("acceptor_resname",  &HBond::acceptor_resname)
        .def_readonly("acceptor_chain",    &HBond::acceptor_chain)
        .def_readonly("acceptor_name",     &HBond::acceptor_name)
        .def_readonly("acceptor_xyz",      &HBond::acceptor_xyz)
        .def_readonly("dist_DA",           &HBond::dist_DA)
        .def_readonly("dist_HA",           &HBond::dist_HA)
        .def_readonly("angle_DHA",         &HBond::angle_DHA)
        .def("__repr__", [](const HBond &hb) {
            return "<HBond " +
                hb.donor_chain + ":" + hb.donor_resname + std::to_string(hb.donor_resid) +
                " " + hb.donor_heavy_name + "-H···" +
                hb.acceptor_chain + ":" + hb.acceptor_resname + std::to_string(hb.acceptor_resid) +
                " " + hb.acceptor_name +
                " DA=" + std::to_string(hb.dist_DA).substr(0,5) +
                " HA=" + std::to_string(hb.dist_HA).substr(0,5) +
                " ang=" + std::to_string(hb.angle_DHA).substr(0,5) + ">";
        });

    // compute_hbonds: donors and acceptors are Model objects (sub-selections);
    // full_model is the complete frame needed for backbone H reconstruction.
    m.def("compute_hbonds",
        [](const Model &donor_model,
           const Model &acceptor_model,
           const Model &full_model,
           float dist_DA_cutoff,
           float dist_HA_cutoff,
           float angle_cutoff) {
            return compute_hbonds(donor_model, acceptor_model, full_model,
                                  dist_DA_cutoff, dist_HA_cutoff, angle_cutoff);
        },
        py::arg("donor_model"),
        py::arg("acceptor_model"),
        py::arg("full_model"),
        py::arg("dist_DA_cutoff") = 3.5f,
        py::arg("dist_HA_cutoff") = 2.5f,
        py::arg("angle_cutoff")   = 90.0f,
        R"doc(
Compute hydrogen bonds between two selections using Baker & Hubbard geometric criteria.

Parameters
----------
donor_model : Model
    Model containing potential donor atoms (a subselection of a frame).
acceptor_model : Model
    Model containing potential acceptor atoms (a subselection of a frame).
full_model : Model
    The complete frame used to reconstruct backbone N-H positions.
dist_DA_cutoff : float, optional
    Maximum donor-heavy to acceptor distance in Å (default 3.5).
dist_HA_cutoff : float, optional
    Maximum hydrogen to acceptor distance in Å (default 2.5).
angle_cutoff : float, optional
    Minimum D-H···A angle in degrees (default 90).

Returns
-------
list[HBond]
    List of detected hydrogen bonds with full geometry information.
)doc");

}
