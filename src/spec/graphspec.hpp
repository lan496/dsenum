#ifndef PYZDD_GRAPHSPEC_H
#define PYZDD_GRAPHSPEC_H

#include <memory>
#include <tdzdd/DdSpec.hpp>
#include "../graph.hpp"

namespace graph {
class SimPath: public tdzdd::PodArrayDdSpec<SimPath, bool, 2> {
    /*
    ZDD to represent simple paths

    let u be Vertex,
    mate[u] = the other endpoint if u is a endpoint of a path
    mate[u] = 0 if u is in the middle of a path
    mate[u] = u (otherwise)

    see TAOCP, Exercise 225, sec.7.1.4, vol.4.1B

    */
    std::shared_ptr<GraphAuxiliary> graphaux_ptr;
public:
    SimPath(const GraphAuxiliary& graphaux) {
        graphaux_ptr = std::make_shared<GraphAuxiliary>(graphaux);
        setArraySize(graphaux.number_of_vertices());
    }
};

}

#endif  // PYZDD_GRAPHSPEC_H
