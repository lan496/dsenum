#include <pybind11/pybind11.h>
#include <tdzdd/DdSpec.hpp>
#include <tdzdd/DdStructure.hpp>
#include "spec/universe.hpp"
#include "spec/combination.hpp"

namespace py = pybind11;

template<typename S>
class SpecBase {};

template<typename S>
class Spec: public SpecBase<S> {
public:
    Spec() {}
};

// A -> Spec<A> -> SpecBase<A>
class A: public Spec<A> {
public:
    A() {}
};

class Structure {
public:
    Structure() {}

    template<typename S>
    void func(SpecBase<S> &spec) {
        return;
    }
};


PYBIND11_MODULE(_pyzdd, m) {
    m.doc() = "wrapper to TdZdd";

    // specifications
    py::class_<tdzdd::DdSpecBase<pyzdd::Combination,2>> PyDdSpecBaseCombination(m, "DdSpecBaseCombination");
    py::class_<tdzdd::DdSpec<pyzdd::Combination,int,2>> PyDdSpecCombination(m, "DdSpecCombination", PyDdSpecBaseCombination);
    py::class_<pyzdd::Combination> PyCombination(m, "Combination", PyDdSpecCombination);
    PyCombination.def(py::init<int, int>(),
                      "specification for DD representing k-combinations out of n items",
                      py::arg("n"),
                      py::arg("k"));

    // DD structure
    py::class_<tdzdd::DdStructure<2>>(m, "Universe")
        .def(py::init<int, bool>(),
             "Universe DD construction",
             py::arg("n"), py::arg("useMP") = false)
        .def("zddSubset", &tdzdd::DdStructure<2>::zddSubset<pyzdd::Combination>);
}
