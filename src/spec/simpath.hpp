#ifndef PYZDD_GRAPHSPEC_H
#define PYZDD_GRAPHSPEC_H

#include <algorithm>
#include <memory>
#include <climits>
#include <tdzdd/DdSpec.hpp>
#include "../graph.hpp"
#include "../type.hpp"

namespace pyzdd {
namespace graph {
namespace simpath {

using Mate = int;

enum MateConstant {
    MIDWAY = -1,
    UNUSED = -2,
};

/// ZDD to represent simple paths from s to t
/// if u belongs to a path, mate[u] = the other endpoint of the path,
/// else if u does not belong to a path, mate[u] = u,
/// else if u is in the midway of a path, mate[u] = MateConstant::MIDWAY
class SimPath: public tdzdd::PodArrayDdSpec<SimPath, Mate, 2> {
    const Vertex V_;
    const Vertex s_;
    const Vertex t_;
    const int E_;
    const int max_frontier_size_;
    const GraphAuxiliary graphaux;

    void initialize(Mate* state) const {
        for (int i = 0; i < max_frontier_size_; ++i) {
            state[i] = MateConstant::UNUSED;
        }
    }

    Mate get_mate(Mate* state, Vertex u) const {
        return state[graphaux.map_vertex_to_position(u)];
    }

    void set_mate(Mate* state, Vertex u, Mate val) const {
        state[graphaux.map_vertex_to_position(u)] = val;
    }

    size_t get_deg(Mate* state, Vertex u) const {
        auto mate_u = get_mate(state, u);
        if (mate_u == u) {
            // isolated
            return 0;
        } else if (mate_u == MateConstant::MIDWAY) {
            // midway of path
            return 2;
        } else {
            // endpoint of a path
            return 1;
        }
    }

    bool is_endpoint(Vertex u) const {
        return ((u == s_) || (u == t_));
    }

public:
    SimPath(const GraphAuxiliary& graphaux, Vertex s, Vertex t) :
        V_(graphaux.number_of_vertices()),
        s_(s),
        t_(t),
        E_(graphaux.number_of_edges()),
        max_frontier_size_(graphaux.get_max_frontier_size()),
        graphaux(graphaux)
    {
        assert(std::is_pod<Mate>::value);
        if (graphaux.number_of_vertices() > INT_MAX) {
            std::cerr << "The number of vertices should be smaller than " << INT_MAX << std::endl;
            exit(1);
        }
        if (s == t) {
            std::cerr << "The endpoints s and t should be different." << std::endl;
            exit(1);
        }

        setArraySize(graphaux.get_max_frontier_size());
    }

    int getRoot(Mate* state) const {
        initialize(state);
        return E_;
    }

    int getChild(Mate* state, int level, int value) const {
        assert((level >= 1) && (level <= E_));

        const InternalEdgeId eid = E_ - level;
        const auto e = graphaux.get_edge(eid);

        // initialize Mate for introduced vertices
        const std::vector<Vertex>& introduced = graphaux.get_introduced(eid);
        for (auto u: introduced) {
            set_mate(state, u, u);
        }

        /*
        std::cerr << std::endl;
        std::cerr << "# call eid=" << eid << ", value=" << value << std::endl;
        std::cerr << "before processing edge" << std::endl;
        print_state(state, level);
        */

        // update state
        const std::vector<Vertex>& frontier = graphaux.get_frontier(eid);
        if (value == 1) {  // take the edge
            auto deg_src = get_deg(state, e.src);
            auto deg_dst = get_deg(state, e.dst);
            auto mate_src = get_mate(state, e.src);
            auto mate_dst = get_mate(state, e.dst);

            // branch occurs
            if ((deg_src == 2) || (deg_dst == 2)) {
                return Terminal::REJECT;
            }
            // degrees of s and t should be one
            if (is_endpoint(e.src) && (deg_src == 1)) {
                return Terminal::REJECT;
            }
            if (is_endpoint(e.dst) && (deg_dst == 1)) {
                return Terminal::REJECT;
            }
            // cycle occurs
            if ((mate_src == e.dst) && (mate_dst == e.src)) {
                return Terminal::REJECT;
            }
            if (((mate_src == s_) && (mate_dst == t_)) || ((mate_src == t_) && (mate_dst == s_))) {
                // s ~~~ e.src -e- e.dst ~~~ t or t ~~~ e.src -e- e.dst ~~~ s
                for (auto u : frontier) {
                    if (!is_endpoint(u) && (u != e.src) && (u != e.dst)) {
                        // vertex u is dangling
                        if (get_deg(state, u) == 1) {
                            return Terminal::REJECT;
                        }
                    }
                }
                // here s-t path is completed
                return Terminal::ACCEPT;
            }

            // process edge
            // other endpoints may not belong to the frontier
            auto other_endpoint_src = get_mate(state, e.src);
            auto other_endpoint_dst = get_mate(state, e.dst);
            assert(other_endpoint_src != MateConstant::UNUSED);
            assert(other_endpoint_src != MateConstant::MIDWAY);
            assert(other_endpoint_dst != MateConstant::UNUSED);
            assert(other_endpoint_dst != MateConstant::MIDWAY);

            if (get_deg(state, e.src) == 0) {
                set_mate(state, e.src, other_endpoint_dst);
            } else {
                set_mate(state, e.src, MateConstant::MIDWAY);
                for (auto u: frontier) {
                    if (u == other_endpoint_src) {
                        set_mate(state, u, other_endpoint_dst);
                    }
                }
            }

            if (get_deg(state, e.dst) == 0) {
                set_mate(state, e.dst, other_endpoint_src);
            } else {
                set_mate(state, e.dst, MateConstant::MIDWAY);
                for (auto u: frontier) {
                    if (u == other_endpoint_dst) {
                        set_mate(state, u, other_endpoint_src);
                    }
                }
            }
        }

        // forget and branch on determined vertex
        const std::vector<Vertex>& forgotten = graphaux.get_forgotten(eid);
        for (auto u : forgotten) {
            size_t deg_u = get_deg(state, u);
            if (is_endpoint(u)) {
                // degrees of s and t should be one
                if (deg_u != 1) {
                    return Terminal::REJECT;
                }
            } else if (deg_u == 1) {
                // degrees of vertices other than s and t should be 0 or 2
                return Terminal::REJECT;
            }

            // forget
            set_mate(state, u, MateConstant::UNUSED);
        }

        /*
        std::cerr << "after update" << std::endl;
        print_state(state, level);
        */

        // if (level == 1) return Terminal::REJECT;
        return level - 1;
    }

    void print_state(Mate* state, int level) const {
        InternalEdgeId eid = E_ - level;

        std::cerr << "frontier" << std::endl;
        const auto& frontier = graphaux.get_frontier(eid);
        for (auto u: frontier) {
            std::cerr << " " << u;
        }
        std::cerr << std::endl;

        std::cerr << "mate" << std::endl;
        for (auto u: frontier) {
            std::cerr << " " << get_mate(state, u);
        }
        std::cerr << std::endl;
    }
};

} // namespace simpath
} // namespace graph
} // namespace pyzdd

#endif  // PYZDD_GRAPHSPEC_H
