#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "Model.h"
#include "Coor.h"

namespace py = pybind11;

PYBIND11_MODULE(core, m) {
    py::class_<Model>(m, "Model")
        .def(py::init<>())
        .def("addAtom", &Model::addAtom)
        .def("clear", &Model::clear)
        .def("size", &Model::size)
        // Add more methods as needed
        ;
    py::class_<Coor>(m, "Coor")
        .def(py::init<>())
        .def("read", &Coor::read)
        .def("write", &Coor::write)
        .def("clear", &Coor::clear)
        .def("add_Model", &Coor::add_Model)
        .def("get_models", &Coor::get_models)
        // Add more methods as needed
        ;
}