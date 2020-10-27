#ifndef PYZDD_GRAPHSPEC_H
#define PYZDD_GRAPHSPEC_H

#include <memory>
#include <tdzdd/DdSpec.hpp>
#include "../graph.hpp"
#include "../type.hpp"

namespace graph {

class FrontierData {
public:
    // degree of the vertex
    Vertex deg;
    // id of connected components contains the vertex
    Vertex comp;
};

// ZDD to represent simple paths from s to t
class SimpleSTPath: public tdzdd::PodArrayDdSpec<SimpleSTPath, FrontierData, 2> {
    Vertex V_;
    Vertex s_;
    Vertex t_;
    int E_;
    int max_frontier_size_;
    std::shared_ptr<GraphAuxiliary> graphaux_ptr;

    void initialize(FrontierData* state) const {
        for (int i = 0; i < max_frontier_size_; ++i) {
            state[i].deg = V_;
            state[i].comp = V_;
        }
    }

    Vertex get_deg(FrontierData* state, Vertex u) const {
        return state[graphaux_ptr->map_vertex_to_position(u)].deg;
    }

    void set_deg(FrontierData* state, Vertex u, Vertex deg) const {
        state[graphaux_ptr->map_vertex_to_position(u)].deg = deg;
    }

    Vertex get_comp(FrontierData* state, Vertex u) const {
        return state[graphaux_ptr->map_vertex_to_position(u)].comp;
    }

    void set_comp(FrontierData* state, Vertex u, Vertex comp) const {
        state[graphaux_ptr->map_vertex_to_position(u)].comp = comp;
    }

    bool is_invalid_degree(FrontierData* state, Vertex u) const {
        if ((u == s_) || (u == t_)) {
            if (get_deg(state, u) > 1) {
                // degrees of s and t should be one
                return true;
            }
        } else if (get_deg(state, u) > 2) {
            // degrees of vertices other than s and t should be 0 or 2
            return true;
        }
        return false;
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
        if (s == t) {
            std::cerr << "The endpoints s and t should be different." << std::endl;
            exit(1);
        }

        setArraySize(graphaux_ptr->get_max_frontier_size());
    }

    int getRoot(FrontierData* state) const {
        initialize(state);
        return E_;
    }

    int getChild(FrontierData* state, int level, int value) const {
        InternalEdgeId eid = E_ - level;
        const Edge& e = graphaux_ptr->get_edge(eid);

        // initialize FrontierData for introduced vertices
        const std::vector<Vertex>& introduced = graphaux_ptr->get_introduced(eid);
        for (const auto u: introduced) {
            set_deg(state, u, 0);
            // initially introduced vertex is an isolated component
            set_comp(state, u, u);
        }

        std::cerr << std::endl;
        std::cerr << "# call eid=" << eid << ", value=" << value << std::endl;
        std::cerr << "before processing edge" << std::endl;
        print_state(state, level);

        // update state
        const std::vector<Vertex>& frontier = graphaux_ptr->get_frontier(eid);
        if (value == 1) {  // take the edge
            set_deg(state, e.src, get_deg(state, e.src) + 1);
            set_deg(state, e.dst, get_deg(state, e.dst) + 1);

            Vertex cs = get_comp(state, e.src);
            Vertex cd = get_comp(state, e.dst);
            if (cs != cd) {
                // merge components cs and cd
                Vertex cmin = std::min(cs, cd);
                Vertex cmax = std::max(cs, cd);
                for (auto u: frontier) {
                    if (get_comp(state, u) == cmax) {
                        set_comp(state, u, cmin);
                    }
                }
            } else {
                // has cycle
                return tdzdd::Terminal::REJECT;
            }
        }

        // branch
        if (is_invalid_degree(state, e.src)) {
            return tdzdd::Terminal::REJECT;
        }
        if (is_invalid_degree(state, e.dst)) {
            return tdzdd::Terminal::REJECT;
        }

        // print_state(state, level, value);
        // forget and branch on determined vertex
        const std::vector<Vertex>& forgotten = graphaux_ptr->get_forgotten(eid);
        for (auto u: forgotten) {
            if ((u == s_) || (u == t_)) {
                // degrees of s and t should be one
                if (get_deg(state, u) != 1) {
                    return tdzdd::Terminal::REJECT;
                }
            } else if ((get_deg(state, u) != 0) && (get_deg(state, u) != 2)) {
                // degrees of vertices other than s and t should be 0 or 2
                return tdzdd::Terminal::REJECT;
            }

            // forget
            set_deg(state, u, V_);
            set_comp(state, u, V_);
        }

        std::cerr << "after update" << std::endl;
        print_state(state, level);

        if (level == 1) return tdzdd::Terminal::ACCEPT;
        return level - 1;
    }

    bool equalTo(const FrontierData* state1, const FrontierData* state2) const {
        return ((state1->deg == state2->deg) && (state1->comp == state2->comp));
    }

    void print_state(FrontierData* state, tdzdd::Level level) const {
        InternalEdgeId eid = E_ - level;

        std::cerr << "[frontier, deg, comp]" << std::endl;
        const auto& frontier = graphaux_ptr->get_frontier(eid);
        for (auto u: frontier) {
            std::cerr << u << " " << get_deg(state, u) << " " << get_comp(state, u) << std::endl;
        }

        std::cerr << "deg :";
        for (FrontierPosition i = 0; i < max_frontier_size_; ++i) {
            std::cerr << " " << state[i].deg;
        }
        std::cerr << std::endl;
        std::cerr << "comp:";
        for (FrontierPosition i = 0; i < max_frontier_size_; ++i) {
            std::cerr << " " << state[i].comp;
        }
        std::cerr << std::endl;
    }
};

}

#endif  // PYZDD_GRAPHSPEC_H
