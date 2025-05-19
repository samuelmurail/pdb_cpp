#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "Model.h"

namespace py = pybind11;

PYBIND11_MODULE(pdb_cpp, m) {
    py::class_<Model>(m, "Model")
        .def(py::init<>())
        .def("read", &Model::read)
        .def("write", &Model::write)
        .def("addAtom", &Model::addAtom)
        .def("clear", &Model::clear)
        .def("size", &Model::size)
        // Add more methods as needed
        ;
}