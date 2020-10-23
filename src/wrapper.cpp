#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <tdzdd/DdSpec.hpp>
#include <tdzdd/DdStructure.hpp>
#include "spec/universe.hpp"
#include "spec/combination.hpp"

namespace py = pybind11;

PYBIND11_MODULE(_pyzdd, m) {
    m.doc() = "wrapper to TdZdd";

    // DD structure
    py::class_<tdzdd::DdStructure<2>> PyDdStructure2(m, "Universe");
    PyDdStructure2.def(py::init<int, bool>(),
                       "Universe DD construction",
                       py::arg("n"), py::arg("useMP") = false);
    PyDdStructure2.def("zddReduce", &tdzdd::DdStructure<2>::zddReduce);
    PyDdStructure2.def("size", &tdzdd::DdStructure<2>::size,
                       "get the number of non-terminal nodes");

    // Set iterator
    using const_iterator = tdzdd::DdStructure<2>::const_iterator;
    py::class_<const_iterator> (m, "const_iterator")
        .def("itemset", &const_iterator::operator*)
        .def("next", &const_iterator::operator++)
        .def(py::self != py::self);
    PyDdStructure2.def("begin", &tdzdd::DdStructure<2>::begin);
    PyDdStructure2.def("end", &tdzdd::DdStructure<2>::end);

    // Combination spec
    py::class_<tdzdd::DdSpecBase<pyzdd::Combination,2>> PyDdSpecBaseCombination(m, "DdSpecBaseCombination");
    py::class_<tdzdd::DdSpec<pyzdd::Combination,int,2>> PyDdSpecCombination(m, "DdSpecCombination", PyDdSpecBaseCombination);
    py::class_<pyzdd::Combination> PyCombination(m, "Combination", PyDdSpecCombination);
    PyCombination.def(py::init<int, int>(),
                      "specification for DD representing k-combinations out of n items",
                      py::arg("n"),
                      py::arg("k"));
    PyDdStructure2.def("zddSubset", &tdzdd::DdStructure<2>::zddSubset<pyzdd::Combination>);
}
