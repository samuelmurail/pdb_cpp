#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include "Model.h"
#include "Coor.h"
#include "TMalign_iface.h"
#include "align.h"
#include "geom.h"
#include "seq_align.h"

namespace py = pybind11;

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
        .def("select_atoms", &Coor::select_atoms, 
            py::arg("selection"), py::arg("frame") = 0, // Specify default value for `frame`
            "Select atoms based on a selection string and an optional frame index")
        .def("get_index_select", &Coor::get_index_select,
            py::arg("selection"), py::arg("frame") = 0, // Specify default value for `frame`
            "Get indices of atoms based on a selection string and an optional frame index")
        .def("select_bool_index", &Coor::select_bool_index)
        .def("get_uniq_chain", &Coor::get_uniq_chain)
        .def("get_aa_sequences", &Coor::get_aa_sequences, 
            py::arg("gap_in_seq") = true, py::arg("frame") = 0, // Specify default values for `gap_in_seq` and `frame`
            "Get the amino acid sequence, optionally including gaps and specifying a frame index")
        // Add more methods as needed
        ;
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
        &sequence_align,  // Directly bind the align function
        py::arg("seq1"),
        py::arg("seq2"),
        py::arg("matrix_file") = "src/pdb_cpp/data/blosum62.txt",
        py::arg("GAP_COST") = -11,
        py::arg("GAP_EXT") = -1,
        py::return_value_policy::take_ownership,  // Python takes ownership of the returned pointer
        "Align two sequences using a substitution matrix and gap penalties");

    // Bind the get_common_atoms function
    m.def("get_common_atoms",
        &get_common_atoms,
        py::arg("coor_1"),
        py::arg("coor_2"),
        py::arg("chain_1") = std::vector<std::string>{"A"},
        py::arg("chain_2") = std::vector<std::string>{"A"},
        py::arg("back_names") = std::vector<std::string>{"C", "N", "O", "CA"},
        py::arg("matrix_file") = "src/pdb_cpp/data/blosum62.txt",
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
        &align_seq_based,
        py::arg("coor_1"),
        py::arg("coor_2"),
        py::arg("chain_1") = std::vector<std::string>{"A"},
        py::arg("chain_2") = std::vector<std::string>{"A"},
        py::arg("back_names") = std::vector<std::string>{"C", "N", "O", "CA"},
        py::arg("matrix_file") = "src/pdb_cpp/data/blosum62.txt",
        py::arg("frame_ref") = 0,
        "Align two coordinate structures using sequence based alignement");
/*
    // Bind the Quaternion struct
    py::class_<Quaternion>(m, "Quaternion")
        .def(py::init<>())
        .def(py::init<float, float, float, float>())
        .def_readwrite("w", &Quaternion::w)
        .def_readwrite("x", &Quaternion::x)
        .def_readwrite("y", &Quaternion::y)
        .def_readwrite("z", &Quaternion::z)
        .def("normalize", &Quaternion::normalize)
        .def("conjugate", &Quaternion::conjugate)
        .def("__mul__", &Quaternion::operator*);

    // Bind quaternion functions
    m.def("makeQ", &makeQ, py::arg("angle"), py::arg("ax"), py::arg("ay"), py::arg("az"),
          "Create a quaternion from axis-angle representation");
    
    m.def("makeW", &makeW, py::arg("x"), py::arg("y"), py::arg("z"),
          "Create a quaternion from x,y,z components, calculating w");
    
    m.def("quaternion_transform", &quaternion_transform, py::arg("q"), py::arg("x"), py::arg("y"), py::arg("z"),
          "Transform a 3D point using quaternion rotation");
    
    m.def("quaternion_rotate", &quaternion_rotate, 
          py::arg("X"), py::arg("Y"),
          "Calculate optimal rotation matrix between two sets of 3D points using quaternion method");
    
    m.def("quaternion_to_matrix", &quaternion_to_matrix, py::arg("q"),
          "Convert quaternion to 3x3 rotation matrix");
    
    m.def("matrix_to_quaternion", &matrix_to_quaternion, py::arg("matrix"),
          "Convert 3x3 rotation matrix to quaternion");
          */
}

// // Alignment align(
// const string &seq1
//  const string &seq2, 
//  onst string &matrix_file,
//  int GAP_COST,
//  int GAP_EXT) {
