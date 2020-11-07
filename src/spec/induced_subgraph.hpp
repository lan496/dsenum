#ifndef PYZDD_INDUCED_SUBGRAPH_H
#define PYZDD_INDUCED_SUBGRAPH_H

#include <type_traits>

#include <tdzdd/DdSpec.hpp>
#include "../graph.hpp"
#include "../type.hpp"

namespace pyzdd {
namespace graph {
namespace induced_subgrah {

/// @brief ZDD to represent vertex-induced subgraphs with the fixed weight.
///        here a weight of a subgraph is defined as the sum of weight of
///        edges which the subgraph possesses.
class VertexInducedSubgraphSpec :
    public tdzdd::HybridDdSpec<VertexInducedSubgraphSpec, Weight, bool, 2> {
    VertexGraphFrontierManager vgfm_;
    Weight target_;
    int max_frontier_size_;
    // the number of vertices
    int V_;
    // weight_sum_lower_bound[vid] is minimum weight sum by taking vid and afterwards
    std::vector<Weight> weight_sum_lower_bound_;
    // weight_sum_upper_bound[vid] is maximum weight sum by taking vid and afterwards
    std::vector<Weight> weight_sum_upper_bound_;
public:
    VertexInducedSubgraphSpec() = delete;
    VertexInducedSubgraphSpec(const VertexInducedSubgraphSpec&) = default;

    VertexInducedSubgraphSpec(const VertexGraphFrontierManager& vgfm, Weight target) :
        vgfm_(vgfm),
        target_(target),
        max_frontier_size_(vgfm.get_max_frontier_size()),
        V_(vgfm.number_of_vertices()) {
        // sanity check on types
        assert(std::is_pod<Weight>::value);
        assert(std::is_pod<bool>::value);

        // max size of `state`
        setArraySize(max_frontier_size_);

        weight_sum_lower_bound_.resize(V_ + 1);
        weight_sum_lower_bound_[V_] = 0;
        for (int vid = V_ - 1; vid >= 0; --vid) {
            weight_sum_lower_bound_[vid] = weight_sum_lower_bound_[vid + 1];
            const std::vector<Edge>& processed_edges = vgfm_.get_processed_edges(vid);
            for (auto e: processed_edges) {
                weight_sum_lower_bound_[vid] += std::min(0, e.weight);
            }
        }

        weight_sum_upper_bound_.resize(V_ + 1);
        weight_sum_upper_bound_[V_] = 0;
        for (int vid = V_ - 1; vid >= 0; --vid) {
            weight_sum_upper_bound_[vid] = weight_sum_upper_bound_[vid + 1];
            const std::vector<Edge>& processed_edges = vgfm_.get_processed_edges(vid);
            for (auto e: processed_edges) {
                weight_sum_upper_bound_[vid] += std::max(0, e.weight);
            }
        }

#ifdef _DEBUG
        std::cerr << "weight lower bound:";
        for (InternalVertexId vid = 0; vid < V_; ++vid) {
            std::cerr << " " << weight_sum_lower_bound_[vid];
        }
        std::cerr << std::endl;
        std::cerr << "weight upper bound:";
        for (InternalVertexId vid = 0; vid < V_; ++vid) {
            std::cerr << " " << weight_sum_upper_bound_[vid];
        }
        std::cerr << std::endl;
#endif
    }

    int getRoot(Weight& current_sum, bool* state) const {
        init_weight_sum(current_sum);
        init_state(state);
        return V_;
    }

    /// `state` remembers selected vertices
    int getChild(Weight& current_sum, bool* state, int level, int value) const {
        InternalVertexId vid = V_ - level;
        Vertex v = vgfm_.get_vertex(vid);

        if (current_sum + get_lower_bound_weight(vid) > target_) {
            // impossible to make the sum to be less than or equal to target
            return Terminal::REJECT;
        }
        if (current_sum + get_upper_bound_weight(vid) < target_) {
            // impossible to make the sum to be greater than or equal to target
            return Terminal::REJECT;
        }

        // initialize state for introduced vertex
        reset_state(state, v);

#ifdef _DEBUG
        std::cerr << std::endl;
        std::cerr << "# vid=" << vid << ", vertex=" << v << ", value=" << value << std::endl;
        std::cerr << "Before processing vertex" << std::endl;
        dump(std::cerr, current_sum, state, level);
#endif

        // iff value == true, choose v in subgraph
        set_state(state, v, value);
        if (value == 1) {
            const std::vector<Edge>& proccessed_edges = vgfm_.get_processed_edges(vid);
            for (Edge e: proccessed_edges) {
                assert(get_state(state, e.src));
                if (get_state(state, e.dst)) {
                    current_sum += e.weight;
                }
            }
        }

#ifdef _DEBUG
        std::cerr << "After processing edges" << std::endl;
        dump(std::cerr, current_sum, state, level);
#endif

        // forget
        const std::vector<Vertex>& forgotten = vgfm_.get_forgotten(vid);
        for (auto uu: forgotten) {
            // vertex uu is no more needed
            reset_state(state, uu);
        }

#ifdef _DEBUG
        std::cerr << "After forgetting vertices" << std::endl;
        dump(std::cerr, current_sum, state, level);
#endif

        if (level == 1) {
            if (current_sum == target_) {
                return Terminal::ACCEPT;
            } else {
                return Terminal::REJECT;
            }
        } else {
            return level - 1;
        }
    }

private:
    void init_weight_sum(Weight& current_sum) const {
        current_sum = 0;
    }

    void init_state(bool* state) const {
        for (int i = 0; i < max_frontier_size_; ++i) {
            state[i] = false;
        }
    }

    bool get_state(bool* state, Vertex v) const {
        return state[vgfm_.map_vertex_to_position(v)];
    }

    void set_state(bool* state, Vertex v, bool value) const {
        state[vgfm_.map_vertex_to_position(v)] = value;
    }

    void reset_state(bool* state, Vertex v) const {
        state[vgfm_.map_vertex_to_position(v)] = false;
    }

    /// @brief minimum attainable sum just before processing the vid-th vertex
    Weight get_lower_bound_weight(InternalVertexId vid) const {
        return weight_sum_lower_bound_[vid];
    }

    /// @brief maximum attainable sum just before processing the vid-th vertex
    Weight get_upper_bound_weight(InternalVertexId vid) const {
        return weight_sum_upper_bound_[vid];
    }

    void dump(std::ostream& os, Weight& current_sum, bool* state, int level) const {
        InternalVertexId vid = V_ - level;
        os << "     weight sum:" << current_sum << std::endl;

        os << "     frontier:";
        const std::vector<Vertex>& frontier = vgfm_.get_frontier(vid);
        for (Vertex uu: frontier) {
            os << " " << uu;
        }
        os << std::endl;

        os << "     state:";
        for (Vertex uu: frontier) {
            os << " " << get_state(state, uu);
        }
        os << std::endl;
    }

};

} // namespace induced_subgraph
} // namespace graph
} // namespace pyzdd

#endif // PYZDD_INDUCED_SUBGRAPH_H
