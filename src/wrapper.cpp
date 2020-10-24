#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <tdzdd/DdSpec.hpp>
#include <tdzdd/DdStructure.hpp>

#include "iterator.hpp"
#include "spec/combination.hpp"
#include "spec/choice.hpp"

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
        .def("next", &const_iterator::operator++)
        .def(py::self != py::self);
    PyDdStructure2.def("begin", &tdzdd::DdStructure<2>::begin);
    PyDdStructure2.def("end", &tdzdd::DdStructure<2>::end);
    m.def("variable_choice", &tdzdd::variable_choice);

    // Specifications
    // Combination spec
    py::class_<tdzdd::DdSpecBase<tdzdd::Combination,2>> PyDdSpecBaseCombination(m, "DdSpecBaseCombination");
    py::class_<tdzdd::DdSpec<tdzdd::Combination,int,2>> PyDdSpecCombination(m, "DdSpecCombination", PyDdSpecBaseCombination);
    py::class_<tdzdd::Combination> PyCombination(m, "Combination", PyDdSpecCombination);
    PyCombination.def(py::init<int, int>(),
                      R"doc(
                        specification for DD representing k-combinations out of n items

                        Parameters
                        ----------
                        n: int
                            the nubmer of items
                        k: int
                            the number of selected items
                      )doc",
                      py::arg("n"),
                      py::arg("k"));
    PyDdStructure2.def("zddSubset", &tdzdd::DdStructure<2>::zddSubset<tdzdd::Combination>);

    // Choice spec
    py::class_<tdzdd::DdSpecBase<tdzdd::Choice,2>> PyDdSpecBaseChoice(m, "DdSpecBaseChoice");
    py::class_<tdzdd::DdSpec<tdzdd::Choice,int,2>> PyDdSpecChoice(m, "DdSpecChoice", PyDdSpecBaseChoice);
    py::class_<tdzdd::Choice> PyChoice(m, "Choice", PyDdSpecChoice);
    PyChoice.def(py::init<int, int, std::vector<int>&>(),
                 R"doc(
                    specification for DD representing combinations out of `n` items such that the number of selected ones in `v` is `k`

                    Parameters
                    ----------
                    n: int
                        the nubmer of items
                    k: int
                        the number of selected items
                    v: List[int]
                        List of 0-indexed variables
                        For example, Choice(4, 2, [0, 2, 3]) represents {{0, 2}, {0, 3}, {2, 3}, {0, 1, 2}, {0, 1, 3}, {1, 2, 3}}.
                 )doc",
                 py::arg("n"),
                 py::arg("k"),
                 py::arg("v"));
    PyDdStructure2.def("zddSubset", &tdzdd::DdStructure<2>::zddSubset<tdzdd::Choice>);
}
