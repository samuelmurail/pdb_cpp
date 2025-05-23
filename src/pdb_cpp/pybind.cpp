#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include "Model.h"
#include "Coor.h"

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
        .def(py::init<>())
        .def(py::init<const std::string&>())  // constructor from filename
        .def("read", &Coor::read)
        .def("write", &Coor::write)
        .def("clear", &Coor::clear)
        .def("size", &Coor::size)
        .def("add_Model", &Coor::add_Model)
        .def("get_Models", &Coor::get_Models)
        .def("select_atoms", &Coor::select_atoms, 
            py::arg("selection"), py::arg("frame") = 0, // Specify default value for `frame`
            "Select atoms based on a selection string and an optional frame index")
        .def("select_bool_index", &Coor::select_bool_index)
        // Add more methods as needed
        ;

}