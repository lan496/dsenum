#ifndef PYZDD_SPANNING_TREE_H
#define PYZDD_SPANNING_TREE_H

#include <algorithm>
#include <memory>
#include <climits>
#include <tdzdd/DdSpec.hpp>
#include "../graph.hpp"
#include "../type.hpp"

namespace pyzdd {
namespace graph {
namespace spanning_forest {

class FrontierData {
public:
    // degree of the vertex
    Vertex deg;
    // id of connected components contains the vertex
    Vertex comp;
};

const Vertex UNUSED = -1;

/// ZDD to represent spanning trees
/// state[u] indicates the index of a connected component that u belongs to.
class SpanningForestSpec: public tdzdd::PodArrayDdSpec<SpanningForestSpec, FrontierData, 2> {
    const Vertex V_;
    const int E_;
    const int max_frontier_size_;
    const GraphAuxiliary graphaux;
public:
    /// @param graphaux frontier manager of a graph
    SpanningForestSpec(const GraphAuxiliary& graphaux) :
        V_(graphaux.number_of_vertices()),
        E_(graphaux.number_of_edges()),
        max_frontier_size_(graphaux.get_max_frontier_size()),
        graphaux(graphaux)
    {
        assert(std::is_pod<FrontierData>::value);
        if (graphaux.number_of_vertices() > std::numeric_limits<Vertex>::max()) {
            std::cerr << "The number of vertices should be smaller than "
                      << std::numeric_limits<Vertex>::max() << std::endl;
            exit(1);
        }

        setArraySize(graphaux.get_max_frontier_size());
    }

    int getRoot(FrontierData* state) const {
        for (int i = 0; i < max_frontier_size_; ++i) {
            state[i].deg = 0;
            state[i].comp = 0;
        }
        return E_;
    }

    int getChild(FrontierData* state, int level, int value) const {
        assert((level >= 1) && (level <= E_));

        const InternalEdgeId eid = E_ - level;
        const auto e = graphaux.get_edge(eid);

        // initialize Mate for introduced vertices
        const std::vector<Vertex>& introduced = graphaux.get_introduced(eid);
        for (auto u: introduced) {
            set_deg(state, u, 0);
            set_comp(state, u, static_cast<Vertex>(u));
        }

#ifdef _DEBUG
        std::cerr << std::endl;
        std::cerr << "# call eid=" << eid << ", value=" << value << std::endl;
        std::cerr << "before processing edge" << std::endl;
        dump_state(std::cerr, state, level);
#endif

        // update state
        const std::vector<Vertex>& frontier = graphaux.get_frontier(eid);
        if (value == 1) {  // take the edge
            auto comp_src = get_comp(state, e.src);
            auto comp_dst = get_comp(state, e.dst);

            // cycle occurs
            if (comp_src == comp_dst) {
                return Terminal::REJECT;
            }

            // process edge
            set_deg(state, e.src, get_deg(state, e.src) + 1);
            set_deg(state, e.dst, get_deg(state, e.dst) + 1);

            // here, comp_src != comp_dst
            auto cmin = std::min(comp_src, comp_dst);
            auto cmax = std::max(comp_src, comp_dst);
            for (auto uf: frontier) {
                if (get_comp(state, uf) == cmin) {
                    // choice cmax so that ComponentId does not decrease.
                    set_comp(state, uf, cmax);
                }
            }
        }

#ifdef _DEBUG
        std::cerr << "after update" << std::endl;
        dump_state(std::cerr, state, level);
#endif

        // branch on determined vertex
        const std::vector<Vertex>& forgotten = graphaux.get_forgotten(eid);
        for (auto v: forgotten) {
            if (get_deg(state, v) == 0) {
                // degree should be at least 1
                return Terminal::REJECT;
            }

            bool comp_found = false;
            for (auto w: frontier) {
                // skip oneself
                if (w == v) {
                    continue;
                }

                // if w has the same component number with v
                if (get_comp(state, w) == get_comp(state, v)) {
                    comp_found = true;
                }
                if (comp_found) {
                    break;
                }
            }
            // if there is no vertex that has the same component with v,
            // the component number of v is determined.
            if (!comp_found) {
                return Terminal::REJECT;
            }
        }

        // forget
        for (auto u: forgotten) {
            set_deg(state, u, 0);
            set_comp(state, u, UNUSED);
        }

#ifdef _DEBUG
        std::cerr << "after forgetting" << std::endl;
        dump_state(std::cerr, state, level);
#endif

        if (level == 1) {
            return Terminal::ACCEPT;
        }
        return level - 1;
    }

    void dump_state(std::ostream& os, FrontierData* state, int level) const {
        InternalEdgeId eid = E_ - level;

        os << "     frontier:";
        const auto& frontier = graphaux.get_frontier(eid);
        for (auto u: frontier) {
            os << " " << u;
        }
        os << std::endl;

        os << "     deg     :";
        for (auto u: frontier) {
            os << " " << get_deg(state, u);
        }
        os << std::endl;

        os << "     comp    :";
        for (auto u: frontier) {
            os << " " << get_comp(state, u);
        }
        os << std::endl;
    }
private:
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
};

} // namespace spanning_tree
} // namespace graph
} // namesapce pyzdd

#endif // PYZDD_SPANNING_TREE_H
