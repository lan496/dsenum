#ifndef PYZDD_GRAPHSPEC_H
#define PYZDD_GRAPHSPEC_H

#include <memory>
#include <tdzdd/DdSpec.hpp>
#include "../graph.hpp"
#include "../type.hpp"

namespace graph {

struct FrontierData {
    // degree of the vertex
    Vertex deg;
    // id of connected components contains the vertex
    Vertex comp;
};

// ZDD to represent simple paths from s to t
class SimpleSTPath: public tdzdd::PodArrayDdSpec<SimpleSTPath, bool, 2> {
    Vertex V_;
    Vertex s_;
    Vertex t_;
    int E_;
    int max_frontier_size_;
    std::shared_ptr<GraphAuxiliary> graphaux_ptr;

    void initialize(FrontierData* state) const {
        for (int i = 0; i < max_frontier_size_; ++i) {
            state[i].deg = 0;
            state[i].comp = 0;
        }
    }

public:
    SimpleSTPath(const GraphAuxiliary& graphaux, Vertex s, Vertex t) :
        V_(graphaux.number_of_vertices()),
        s_(s),
        t_(t),
        E_(graphaux.number_of_edges()),
        max_frontier_size_(graphaux.get_max_frontier_size())
    {
        graphaux_ptr = std::make_shared<GraphAuxiliary>(graphaux);

        if (graphaux_ptr->number_of_vertices() > SHRT_MAX) {
            std::cerr << "The number of vertices should be smaller than " << SHRT_MAX << std::endl;
            exit(1);
        }

        setArraySize(graphaux_ptr->get_max_frontier_size());
    }

    int getRoot(FrontierData* state) const {
        initialize(state);
        return E_;
    }

    int getChild(FrontierData* state, tdzdd::Level level, int value) const {
        InternalEdgeId eid = E_ - level;
        const auto& e = graphaux_ptr->get_edge(eid);

        // initialize FrontierData for introduced vertices
    }
};

}

#endif  // PYZDD_GRAPHSPEC_H
