#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include "Model.h"
#include "Coor.h"
#include "TMAlign.h"
#include "sequence_align.h"

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
        .def("select_bool_index", &Coor::select_bool_index)
        .def("get_uniq_chain", &Coor::get_uniq_chain)
        .def("get_aa_sequences", &Coor::get_aa_sequences, 
            py::arg("gap_in_seq") = true, py::arg("frame") = 0, // Specify default values for `gap_in_seq` and `frame`
            "Get the amino acid sequence, optionally including gaps and specifying a frame index")
        // Add more methods as needed
        ;
    // Bind the TMAlign function
    m.def("compute_SS",
        static_cast<std::vector<std::vector<std::string>>(*)(const Coor&, bool)>(&compute_SS),
        py::arg("coor"),
        py::arg("gap_in_seq") = false,
        "Compute secondary structure for all models in a Coor object");

    // Bind the Alignment structure
    py::class_<Alignment>(m, "Alignment")
        .def(py::init<>())  // Default constructor
        .def_readwrite("seq1", &Alignment::seq1, "Aligned sequence 1")
        .def_readwrite("seq2", &Alignment::seq2, "Aligned sequence 2")
        .def_readwrite("score", &Alignment::score, "Alignment score");

    // Bind the align function
    m.def("sequence_align",
        &sequence_align,  // Directly bind the align function
        py::arg("seq1"),
        py::arg("seq2"),
        py::arg("matrix_file") = "src/pdb_cpp/data/blosum62.txt",
        py::arg("GAP_COST") = -11,
        py::arg("GAP_EXT") = -1,
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
}

// // Alignment align(
// const string &seq1
//  const string &seq2, 
//  onst string &matrix_file,
//  int GAP_COST,
//  int GAP_EXT) {
