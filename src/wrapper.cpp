#include <pybind11/pybind11.h>
#include <tdzdd/DdSpec.hpp>
#include <tdzdd/DdStructure.hpp>
#include "spec/universe.hpp"
#include "spec/combination.hpp"

namespace py = pybind11;

PYBIND11_MODULE(_pyzdd, m) {
    m.doc() = "wrapper to TdZdd";

    // specifications
    py::class_<pyzdd::Combination>(m, "Combination")
        .def(py::init<int, int>(),
             "specification for DD representing k-combinations out of n items",
             py::arg("n"),
             py::arg("k"));

    // DD structure
    py::class_<tdzdd::DdStructure<2>>(m, "Universe")
        .def(py::init<int, bool>(),
             "Universe DD construction",
             py::arg("n"), py::arg("useMP") = false)
        .def("zddSubset", &pyzdd::DdStructure<2>::zddSubset);
}
