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
namespace stpath {

class FrontierData {
public:
    // degree of the vertex
    Vertex deg;
    // id of connected components contains the vertex
    Vertex comp;
};

// ZDD to represent simple paths from s to t
class SimpleSTPath: public tdzdd::PodArrayDdSpec<SimpleSTPath, FrontierData, 2> {
    const Vertex V_;
    const Vertex s_;
    const Vertex t_;
    const int E_;
    const int max_frontier_size_;
    const Level s_introduced_level_;
    const Level t_introduced_level_;
    const GraphAuxiliary graphaux;

    void initialize(FrontierData* state) const {
        for (int i = 0; i < max_frontier_size_; ++i) {
            state[i].deg = 0;
            state[i].comp = 0;
        }
    }

    Vertex get_deg(FrontierData* state, Vertex u) const {
        return state[graphaux.map_vertex_to_position(u)].deg;
    }

    void set_deg(FrontierData* state, Vertex u, Vertex deg) const {
        state[graphaux.map_vertex_to_position(u)].deg = deg;
    }

    Vertex get_comp(FrontierData* state, Vertex u) const {
        return state[graphaux.map_vertex_to_position(u)].comp;
    }

    void set_comp(FrontierData* state, Vertex u, Vertex comp) const {
        state[graphaux.map_vertex_to_position(u)].comp = comp;
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
        max_frontier_size_(graphaux.get_max_frontier_size()),
        s_introduced_level_(graphaux.get_vertex_introduced_level(s)),
        t_introduced_level_(graphaux.get_vertex_introduced_level(t)),
        graphaux(graphaux)
    {
        assert(std::is_pod<FrontierData>::value);
        if (graphaux.number_of_vertices() > SHRT_MAX) {
            std::cerr << "The number of vertices should be smaller than " << SHRT_MAX << std::endl;
            exit(1);
        }
        if (s == t) {
            std::cerr << "The endpoints s and t should be different." << std::endl;
            exit(1);
        }

        setArraySize(graphaux.get_max_frontier_size());
    }

    int getRoot(FrontierData* state) const {
        initialize(state);
        return E_;
    }

    int getChild(FrontierData* state, int level, int value) const {
        assert((level >= 1) && (level <= E_));

        InternalEdgeId eid = E_ - level;
        const Edge& e = graphaux.get_edge(eid);

        // initialize FrontierData for introduced vertices
        const std::vector<Vertex>& introduced = graphaux.get_introduced(eid);
        for (const auto v: introduced) {
            set_deg(state, v, 0);
            // initially introduced vertex is an isolated component
            set_comp(state, v, v);
        }

        // std::cerr << std::endl;
        // std::cerr << "# call eid=" << eid << ", value=" << value << std::endl;
        // std::cerr << "before processing edge" << std::endl;
        // print_state(state, level);

        // update state
        const std::vector<Vertex>& frontier = graphaux.get_frontier(eid);
        if (value == 1) {  // take the edge
            set_deg(state, e.src, get_deg(state, e.src) + 1);
            set_deg(state, e.dst, get_deg(state, e.dst) + 1);

            Vertex cs = get_comp(state, e.src);
            Vertex cd = get_comp(state, e.dst);
            if (cs != cd) {
                // merge components cs and cd
                Vertex cmin = std::min(cs, cd);
                Vertex cmax = std::max(cs, cd);
                for (auto v: frontier) {
                    // choice cmax so that `comp` does not decrease.
                    if (get_comp(state, v) == cmin) {
                        set_comp(state, v, cmax);
                    }
                }
            } else {
                // has cycle
                return Terminal::REJECT;
            }
        }

        // branch
        if (is_invalid_degree(state, e.src)) {
            return Terminal::REJECT;
        }
        if (is_invalid_degree(state, e.dst)) {
            return Terminal::REJECT;
        }

        // print_state(state, level, value);
        // forget and branch on determined vertex
        const std::vector<Vertex>& forgotten = graphaux.get_forgotten(eid);
        for (size_t i = 0; i < forgotten.size(); ++i) {
            Vertex v = forgotten[i];
            if ((v == s_) || (v == t_)) {
                // degrees of s and t should be one
                if (get_deg(state, v) != 1) {
                    return Terminal::REJECT;
                }
            } else if ((get_deg(state, v) != 0) && (get_deg(state, v) != 2)) {
                // degrees of vertices other than s and t should be 0 or 2
                return Terminal::REJECT;
            }

            // find frontier vertex whose degree is more than 1 and belong
            // the other connected component than u
            bool comp_found = false;
            bool deg_found = false;
            for (Vertex w_frontier: frontier) {
                // skip oneself
                if (w_frontier == v) {
                    continue;
                }
                // skip forgotten vertices
                bool found_leaved = false;
                for (size_t k = 0; k < i; ++k) {
                    if (w_frontier == forgotten[k]) {
                        found_leaved = true;
                        break;
                    }
                }
                if (found_leaved) {
                    continue;
                }

                if (get_comp(state, w_frontier) == get_comp(state, v)) {
                    comp_found = true;
                }
                if (get_deg(state, w_frontier) > 0) {
                    deg_found = true;
                }
                if (deg_found && comp_found) {
                    break;
                }
            }

            // when a component number of u is determined
            if (!comp_found) {
                assert(get_deg(state, v) <= 2);
                // there are multiple components. contradiction.
                if (get_deg(state, v) > 0 && deg_found) {
                    return Terminal::REJECT;
                } else if (get_deg(state, v) > 0) {
                    // s or t cannot belong to the same components of u. contradiction.
                    if (level > s_introduced_level_ || level > t_introduced_level_) {
                        return Terminal::REJECT;
                    } else {
                        // s-t path already completes
                        return Terminal::ACCEPT;
                    }
                }
            }

            // forget
            set_deg(state, v, V_);
            set_comp(state, v, V_);
        }

        // std::cerr << "after update" << std::endl;
        // print_state(state, level);

        if (level == 1) return Terminal::REJECT;
        return level - 1;
    }

    void print_state(FrontierData* state, Level level) const {
        InternalEdgeId eid = E_ - level;

        std::cerr << "[frontier, deg, comp]" << std::endl;
        const auto& frontier = graphaux.get_frontier(eid);
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

} // namespace stpath
} // namesapce graph
} // namespace pyzdd

#endif  // PYZDD_GRAPHSPEC_H
