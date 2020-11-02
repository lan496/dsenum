#include <vector>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <tdzdd/DdSpec.hpp>
#include <tdzdd/DdStructure.hpp>

#include "iterator.hpp"
#include "graph.hpp"
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
    PyDdStructure2.def("cardinality", &tdzdd::DdStructure<2>::zddCardinality,
                       "count the number of sets in the family represented by this ZDD. Returned type is str because cardinality may be too enormous than int64.");

    // Set iterator
    using const_iterator = tdzdd::DdStructure<2>::const_iterator;
    py::class_<const_iterator> (m, "const_iterator")
        .def("next", &const_iterator::operator++)
        .def(py::self != py::self);
    PyDdStructure2.def("begin", &tdzdd::DdStructure<2>::begin);
    PyDdStructure2.def("end", &tdzdd::DdStructure<2>::end);
    m.def("variable_choice", &pyzdd::variable_choice);

    // Graph
    py::class_<pyzdd::graph::Edge> (m, "Edge")
        .def(py::init<int,int>())
        .def(py::init<int,int,int>())
        .def("__repr__", [](const pyzdd::graph::Edge &e) {
            return "<Edge(src=" + std::to_string(e.src) + ", dst=" + std::to_string(e.dst) + ", weight=" + std::to_string(e.weight) + ")>";
        });
    py::class_<pyzdd::graph::Graph> (m, "Graph");
    py::class_<pyzdd::graph::GraphAuxiliary> (m, "GraphAuxiliary")
        .def(py::init<const pyzdd::graph::Graph&>())
        .def_property_readonly("max_frontier_size", &pyzdd::graph::GraphAuxiliary::get_max_frontier_size)
        .def("edge_order", &pyzdd::graph::GraphAuxiliary::get_edge_order)
        .def("frontier", &pyzdd::graph::GraphAuxiliary::get_frontier)
        .def("introduced", &pyzdd::graph::GraphAuxiliary::get_introduced)
        .def("forgotten", &pyzdd::graph::GraphAuxiliary::get_forgotten)
        .def("map_vertex_to_position", &pyzdd::graph::GraphAuxiliary::map_vertex_to_position);

    // Specifications
    // Combination spec
    py::class_<tdzdd::DdSpecBase<pyzdd::combination::Combination,2>> PyDdSpecBaseCombination(m, "DdSpecBaseCombination");
    py::class_<tdzdd::DdSpec<pyzdd::combination::Combination,int,2>> PyDdSpecCombination(m, "DdSpecCombination", PyDdSpecBaseCombination);
    py::class_<pyzdd::combination::Combination> PyCombination(m, "Combination", PyDdSpecCombination);
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
    PyDdStructure2.def("zddSubset", &tdzdd::DdStructure<2>::zddSubset<pyzdd::combination::Combination>);

    // Choice spec
    py::class_<tdzdd::DdSpecBase<pyzdd::choice::Choice,2>> PyDdSpecBaseChoice(m, "DdSpecBaseChoice");
    py::class_<tdzdd::DdSpec<pyzdd::choice::Choice,int,2>> PyDdSpecChoice(m, "DdSpecChoice", PyDdSpecBaseChoice);
    py::class_<pyzdd::choice::Choice> PyChoice(m, "Choice", PyDdSpecChoice);
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
    PyDdStructure2.def("zddSubset", &tdzdd::DdStructure<2>::zddSubset<pyzdd::choice::Choice>);
}
