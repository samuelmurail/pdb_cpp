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
        .def(py::init<const std::string&>())  // constructor from filename
        .def("read", &Coor::read)
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

}
