#ifndef PYZDD_GRAPH_H
#define PYZDD_GRAPH_H

#include <iostream>
#include <vector>
#include <queue>
#include <set>
#include <algorithm>
#include <limits>
#include <cassert>
#include "type.hpp"

namespace pyzdd {
namespace graph {
// =============================================================================
// Edge and Graph
// =============================================================================
using Vertex = int;
using Weight = int;
/// @brief represent edge of graph
struct Edge {
    Vertex src, dst;
    Weight weight;
    Edge() {}
    Edge(Vertex src, Vertex dst) :
        src(src), dst(dst), weight(1) {}
    Edge(Vertex src, Vertex dst, Weight weight) :
        src(src), dst(dst), weight(weight) {}
};

bool operator==(const Edge& lhs, const Edge& rhs) {
    return (lhs.src == rhs.src) && (lhs.dst == rhs.dst) && (lhs.weight == rhs.weight);
}

std::ostream& operator<<(std::ostream& os, const Edge& e) {
    os << "(src=" << e.src << " ,dst=" << e.dst << ", w=" << e.weight << ")";
    return os;
}

using Edges = std::vector<Edge>;
using Graph = std::vector<Edges>;

/// @brief add undirected edge {u, v} with weight w into graph
void add_undirected_edge(Graph &graph, Vertex u, Vertex v, Weight w) {
    graph[u].emplace_back(Edge(u, v, w));
    graph[v].emplace_back(Edge(v, u, w));
}

int get_number_of_edges(const Graph &undirected_graph) {
    int V = undirected_graph.size();
    int E = 0;
    for (Vertex u = 0; u < static_cast<Vertex>(V); ++u) {
        for (auto& e: undirected_graph[u]) {
            if (e.dst > u) {
                ++E;
            }
        }
    }
    return E;
}

/// @brief check if a given has no multi-edges and self-loops.
bool is_simple_graph(const Graph& graph) {
    Vertex V = graph.size();
    for (Vertex u = 0; u < V; ++u) {
        // check no multiedge
        std::set<Vertex> dsts;
        for (auto& e: graph[u]) {
            dsts.insert(e.dst);
        }
        if (dsts.size() != graph[u].size()) {
            return false;
        }

        // check no self-edge
        for (auto& e: graph[u]) {
            if (e.dst == u) {
                return false;
            }
        }
    }
    return true;
}

// Internal types for variables in DD
using InternalVertexId = Vertex;
// Internal types for variables in DD
using InternalEdgeId = int;
// position in PodArray
using FrontierPosition = Vertex;

/// @brief Frontier manager for ZDD represeting edge-induced subgraphs
/// TODO: rename to EdgeGraphFrontierManager
class GraphAuxiliary {
private:
    /// number of vertices
    size_t V_;
    /// number of edges
    size_t E_;
    /// Edge edge_order[i] is processed at the i-th
    /// edge_order: InternalEdgeId -> Edge
    std::vector<Edge> edge_order_;

    /// Just before processing the i-th edge, the states of vertices of frontiers[i] are required.
    /// For each i, frontiers_[i] is sorted by ascending order.
    std::vector<std::vector<Vertex>> frontiers_;
    /// When processing the i-the edge, vertices introduced[i] enter.
    std::vector<std::vector<Vertex>> introduced_;
    /// After processing the i-th edge, vertices forgotten[i] are no more needed.
    std::vector<std::vector<Vertex>> forgotten_;

    /// mapping_vertex_pos_[v] is position of vertex v in frontiers
    std::vector<FrontierPosition> mapping_vertex_pos_;

    /// the maximum size of frontiers
    int max_frontier_size_;
public:
    GraphAuxiliary() {}
    GraphAuxiliary(const Graph& graph) : V_(graph.size())  {
        // sanity check on graph
        assert(is_simple_graph(graph));
        for (Vertex u = 0; u < static_cast<Vertex>(V_); ++u) {
            if (graph[u].empty()) {
                std::cerr << "TODO: handle isolated vertices" << std::endl;
                exit(1);
            }
        }

        if (V_ > SHRT_MAX) {
            std::cerr << "The number of vertices should be smaller than " << SHRT_MAX << std::endl;
            exit(1);
        }

        // number of edges
        E_ = 0;
        for (Vertex u = 0; u < static_cast<Vertex>(V_); ++u) {
            for (auto& e: graph[u]) {
                if (e.dst > u) {
                    ++E_;
                }
            }
        }

        // order edges by BFS
        order_edges(graph);
        // construct frontiers, intorduced, forgotten, max_frontier_size
        construct_frontiers();
    }

    /// @brief prepare frontiers with a given varible order.
    GraphAuxiliary(const Graph& graph, const std::vector<Edge>& edge_order) :
        V_(graph.size()),
        edge_order_(edge_order)
    {
        assert(is_simple_graph(graph));
        if (V_ > SHRT_MAX) {
            std::cerr << "The number of vertices should be smaller than " << SHRT_MAX << std::endl;
            exit(1);
        }

        // number of edges
        E_ = 0;
        for (Vertex u = 0; u < static_cast<Vertex>(V_); ++u) {
            for (auto& e: graph[u]) {
                if (e.dst > u) {
                    ++E_;
                }
            }
        }

        assert(edge_order.size() == E_);

        // construct frontiers, intorduced, forgotten, max_frontier_size
        construct_frontiers();
    }

    /// @brief check if a given has no multi-edges and self-loops.
    bool is_simple_graph(const Graph& graph) const {
        Vertex V = graph.size();
        for (Vertex u = 0; u < V; ++u) {
            // check no multiedge
            std::set<Vertex> dsts;
            for (auto& e: graph[u]) {
                dsts.insert(e.dst);
            }
            if (dsts.size() != graph[u].size()) {
                return false;
            }

            // check no self-edge
            for (auto& e: graph[u]) {
                if (e.dst == u) {
                    return false;
                }
            }
        }
        return true;
    }

    size_t number_of_vertices() const {
        return V_;
    }

    size_t number_of_edges() const {
        return E_;
    }

    size_t get_max_frontier_size() const {
        return max_frontier_size_;
    }

    const std::vector<Edge>& get_edge_order() const {
        return edge_order_;
    }

    const Edge& get_edge(InternalEdgeId eid) const {
        return edge_order_[eid];
    }

    const std::vector<Vertex>& get_frontier(InternalEdgeId eid) const {
        return frontiers_[eid];
    }

    const std::vector<Vertex>& get_introduced(InternalEdgeId eid) const {
        return introduced_[eid];
    }

    const std::vector<Vertex>& get_forgotten(InternalEdgeId eid) const {
        return forgotten_[eid];
    }

    const FrontierPosition map_vertex_to_position(Vertex u) const {
        return mapping_vertex_pos_[u];
    }

    pyzdd::Level get_vertex_introduced_level(Vertex u) const {
        InternalEdgeId E = number_of_edges();
        for (InternalEdgeId eid = 0; eid < E; ++eid) {
            for (const Vertex &ve: get_introduced(eid)) {
                if (ve == u) return E - eid;
            }
        }
        return 0;  // isolated vertex
    }

    void dump(std::ostream& os) const {
        int V = number_of_vertices();
        int E = number_of_edges();
        os << "V=" << V << std::endl;
        os << "E=" << E << std::endl;

        os << "frontiers" << std::endl;
        for (InternalEdgeId eid = 0; eid < E; ++eid) {
            os << "     " << eid << ":";
            for (const auto& u: get_frontier(eid)) {
                os << " " << u;
            }
            os << std::endl;
        }
        os << "introduced" << std::endl;
        for (InternalEdgeId eid = 0; eid < E; ++eid) {
            os << "     " << eid << ":";
            for (const auto& u: get_introduced(eid)) {
                os << " " << u;
            }
            os << std::endl;
        }
        os << "forgotten" << std::endl;
        for (InternalEdgeId eid = 0; eid < E; ++eid) {
            os << "     " << eid << ":";
            for (const auto& u: get_forgotten(eid)) {
                os << " " << u;
            }
            os << std::endl;
        }
        os << "mapping" << std::endl;
        for (Vertex u = 0; u < static_cast<Vertex>(V_); ++u) {
            os << " " << map_vertex_to_position(u);
        }
        os << std::endl;
    }

    /// @deprecated
    void print() const {
        dump(std::cerr);
    }
private:
    /// @brief determine variable order by BFS
    /// @return bool true iff the graph is connected.
    void order_edges(const Graph& graph) {
        Vertex V = number_of_vertices();
        size_t E = number_of_edges();

        // order edges by BFS
        edge_order_.resize(E);
        InternalEdgeId edge_count = 0;
        std::vector<bool> visited_vertex(V, false);
        for (Vertex u = 0; u < V; ++u) {
            if (visited_vertex[u]) continue;

            std::queue<Vertex> que;
            que.push(u);
            while (!que.empty()) {
                Vertex v = que.front(); que.pop();
                if (visited_vertex[v]) {
                    continue;
                }
                visited_vertex[v] = true;
                for (auto e: graph[v]) {
                    if (!visited_vertex[e.dst]) {
                        edge_order_[edge_count++] = e;
                        que.push(e.dst);
                    }
                }
            }
        }
        assert(edge_count == static_cast<InternalEdgeId>(E));
    }

    /// @brief construct `introduced_` and `forgotten_`
    void construct_introduced_forgotten() {
        Vertex V = number_of_vertices();
        InternalEdgeId E = number_of_edges();
        const auto& edge_order = get_edge_order();

        // calculate introduced
        introduced_.resize(E);
        std::vector<bool> visited_vertex(V, false);
        for (InternalEdgeId eid = 0; eid < E; ++eid) {
            auto e = edge_order[eid];
            if (!visited_vertex[e.src]) {
                introduced_[eid].emplace_back(e.src);
                visited_vertex[e.src] = true;
            }
            if (!visited_vertex[e.dst]) {
                introduced_[eid].emplace_back(e.dst);
                visited_vertex[e.dst] = true;
            }
        }
        for (Vertex u = 0; u < V; ++u) {
            assert(visited_vertex[u]);
        }

        // calculate forgotten by backward
        forgotten_.resize(E);
        std::vector<bool> processed_vertex(V, false);
        for (int eid = E - 1; eid >= 0; eid--) {  // Be careful overflow!
            auto e = edge_order[eid];
            if (!processed_vertex[e.src]) {
                forgotten_[eid].emplace_back(e.src);
                processed_vertex[e.src] = true;
            }
            if (!processed_vertex[e.dst]) {
                forgotten_[eid].emplace_back(e.dst);
                processed_vertex[e.dst] = true;
            }
        }
        for (Vertex u = 0; u < V; ++u) {
            assert(processed_vertex[u]);
        }
    }

    void construct_frontiers() {
        Vertex V = number_of_vertices();
        InternalEdgeId E = number_of_edges();

        construct_introduced_forgotten();

        // calculate frontiers
        frontiers_.resize(E);
        std::set<Vertex> bag;
        for (InternalEdgeId eid = 0; eid < E; ++eid) {
            // introduced vertices
            for (auto v: introduced_[eid]) {
                bag.insert(v);
            }
            // determine frontier
            for (auto v: bag) {
                frontiers_[eid].emplace_back(v);
            }
            // forgotten vertices
            for (auto v: forgotten_[eid]) {
                bag.erase(v);
            }
        }
        assert(bag.empty());

        // frontier size
        max_frontier_size_ = 0;
        for (const auto& f: frontiers_) {
            max_frontier_size_ = std::max(static_cast<int>(f.size()), max_frontier_size_);
        }

        // map vertex to position
        // fill unused entry with V
        mapping_vertex_pos_.assign(V, V);
        std::vector<bool> used(max_frontier_size_, false);
        for (InternalEdgeId eid = 0; eid < E; ++eid) {
            // introduce vertex
            for (auto u: introduced_[eid]) {
                bool success = false;
                for (FrontierPosition i = 0; i < max_frontier_size_; ++i) {
                    if (!used[i]) {
                        mapping_vertex_pos_[u] = i;
                        used[i] = true;
                        success = true;
                        break;
                    }
                }
                assert(success);
            }
            // forget vertex
            for (auto u: forgotten_[eid]) {
                FrontierPosition released = mapping_vertex_pos_[u];
                used[released] = false;
            }
        }
        for (Vertex u = 0; u < V; ++u) {
            assert(mapping_vertex_pos_[u] != V);
        }
        for (FrontierPosition i = 0; i < max_frontier_size_; ++i) {
            assert(!used[i]);
        }
    }
};

/// @brief Frontier manager for ZDD represeting vertex-induced subgraphs
class VertexGraphFrontierManager {
private:
    /// number of vertices
    size_t V_;
    /// number of edges
    size_t E_;
    /// undirected simple graph
    Graph graph_;

    /// vertex_order_[i] is processed at the i-th
    /// vertex_order_: InternalVertexId-> Vertex
    std::vector<Vertex> vertex_order_;
    /// mapping_vertex: Vertex -> InternalVertexId
    std::vector<InternalVertexId> mapping_vertex_;

    /// Just before processing the i-th vertex, the states of vertices of frontiers[i] are required.
    /// For each i, frontiers_[i] is sorted by ascending order.
    std::vector<std::vector<Vertex>> frontiers_;
    /// After processing the i-th vertex, processed_[i] can also be processed.
    /// processed_: InternalVertexId -> [Edge]
    std::vector<std::vector<Edge>> processed_;
    /// After processing the i-th vertex, vertices forgotten[i] are no more needed.
    std::vector<std::vector<Vertex>> forgotten_;

    /// mapping_vertex_pos_[v] is position of vertex v in frontiers
    std::vector<FrontierPosition> mapping_vertex_pos_;

    /// the maximum size of frontiers
    int max_frontier_size_;
public:
    VertexGraphFrontierManager() {}
    VertexGraphFrontierManager(const Graph& graph) :
        V_(graph.size()),
        E_(get_number_of_edges(graph)),
        graph_(graph)
    {
        // sanity check on graph
        assert(is_simple_graph(graph_));

        if (V_ > std::numeric_limits<Vertex>::max()) {
            std::cerr << "The number of vertices should be smaller than " << std::numeric_limits<Vertex>::max() << std::endl;
            exit(1);
        }

        // order vertices by BFS
        order_vertices();
        // construct frontiers, introduced, forgotten, max_frontier_size
        construct_frontiers();
    }

    /// @brief prepare frontiers with a given variable order.
    VertexGraphFrontierManager(const Graph& graph, const std::vector<Vertex>& vertex_order) :
        V_(graph.size()),
        E_(get_number_of_edges(graph)),
        graph_(graph),
        vertex_order_(vertex_order)
    {
        assert(is_simple_graph(graph_));
        if (V_ > std::numeric_limits<Vertex>::max()) {
            std::cerr << "The number of vertices should be smaller than " << std::numeric_limits<Vertex>::max() << std::endl;
            exit(1);
        }

        assert(vertex_order_.size() == V_);

        // construct frontiers, introduced, forgotten, max_frontier_size
        construct_frontiers();
    }

    size_t number_of_vertices() const {
        return V_;
    }

    size_t number_of_edges() const {
        return E_;
    }

    size_t get_max_frontier_size() const {
        return max_frontier_size_;
    }

    const std::vector<Vertex>& get_vertex_order() const {
        return vertex_order_;
    }

    const std::vector<InternalVertexId>& get_mapping_vertex() const {
        return mapping_vertex_;
    }

    /// @brief return proccess vertex at the vid-th in DD
    const Vertex& get_vertex(InternalVertexId vid) const {
        return vertex_order_[vid];
    }

    /// @brief return when v is proccessed
    InternalVertexId map_to_internal_vertex_id(Vertex v) const {
        return mapping_vertex_[v];
    }

    const std::vector<Edge>& get_processed_edges(InternalVertexId vid) const {
        return processed_[vid];
    }

    const std::vector<Vertex>& get_frontier(InternalVertexId vid) const {
        return frontiers_[vid];
    }

    const std::vector<Vertex>& get_forgotten(InternalVertexId vid) const {
        return forgotten_[vid];
    }

    const FrontierPosition map_vertex_to_position(Vertex u) const {
        return mapping_vertex_pos_[u];
    }

    void dump(std::ostream& os) const {
        int V = number_of_vertices();
        int E = number_of_edges();
        os << "V=" << V << std::endl;
        os << "E=" << E << std::endl;

        os << "frontiers" << std::endl;
        for (InternalVertexId vid = 0; vid < V; ++vid) {
            os << "     " << vid << "(v=" << get_vertex(vid) << "):";
            for (const auto& u: get_frontier(vid)) {
                os << " " << u;
            }
            os << std::endl;
        }

        os << "processed edges" << std::endl;
        for (InternalVertexId vid = 0; vid < V; ++vid) {
            os << "     " << vid << "(v=" << get_vertex(vid) << "):";
            for (const auto& e: get_processed_edges(vid)) {
                os << " " << e;
            }
            os << std::endl;
        }

        os << "forgotten" << std::endl;
        for (InternalVertexId vid = 0; vid < V; ++vid) {
            os << "     " << vid << "(v=" << get_vertex(vid) << "):";
            for (const auto& u: get_forgotten(vid)) {
                os << " " << u;
            }
            os << std::endl;
        }

        os << "mapping" << std::endl;
        for (Vertex u = 0; u < static_cast<Vertex>(V_); ++u) {
            os << " " << map_vertex_to_position(u);
        }
        os << std::endl;

        os << "max frontier size: " << get_max_frontier_size() << std::endl;
    }

private:
    /// @brief determine variable order by BFS
    void order_vertices() {
        Vertex V = number_of_vertices();

        // order edges by BFS
        vertex_order_.resize(V);
        InternalVertexId vertex_count = 0;
        std::vector<bool> visited_vertex(V, false);
        for (Vertex u = 0; u < V; ++u) {
            if (visited_vertex[u]) continue;

            std::queue<Vertex> que;
            que.push(u);
            while (!que.empty()) {
                Vertex v = que.front(); que.pop();
                if (visited_vertex[v]) {
                    continue;
                }
                visited_vertex[v] = true;
                vertex_order_[vertex_count++] = v;
                for (auto e: graph_[v]) {
                    if (!visited_vertex[e.dst]) {
                        que.push(e.dst);
                    }
                }
            }
        }
        assert(vertex_count == static_cast<InternalVertexId>(V));

        // mapping_vertex_
        mapping_vertex_.resize(V);
        for (InternalVertexId vid = 0; vid < V; ++vid) {
            mapping_vertex_[vertex_order_[vid]] = vid;
        }
    }

    void construct_frontiers() {
        auto V = number_of_vertices();

        // calculate processed_
        processed_.resize(V);
        int count_processed_edges = 0;
        for (InternalVertexId vid = 0; vid < static_cast<InternalVertexId>(V); ++vid) {
            Vertex v = get_vertex(vid);
            for (Edge e: graph_[v]) {
                if (map_to_internal_vertex_id(e.dst) <= vid) {
                    processed_[vid].emplace_back(e);
                    count_processed_edges++;
                }
            }
        }
        assert(count_processed_edges == static_cast<int>(E_));

        // calculate forgotten vertices
        forgotten_.resize(V);
        for (InternalVertexId vid = 0; vid < static_cast<InternalVertexId>(V); ++vid) {
            Vertex v = get_vertex(vid);
            InternalVertexId forget_v = vid;
            for (Edge e: graph_[v]) {
                forget_v = std::max(map_to_internal_vertex_id(e.dst), forget_v);
            }
            forgotten_[forget_v].emplace_back(v);
        }

        // calculate frontiers
        frontiers_.resize(V);
        std::set<Vertex> bag;
        for (InternalVertexId vid = 0; vid < static_cast<InternalVertexId>(V); ++vid) {
            // introduced vertex
            Vertex v = get_vertex(vid);
            bag.insert(v);

            // determine frontier
            for (auto uu: bag) {
                frontiers_[vid].emplace_back(uu);
            }

            // forgotten vertices
            for (auto uu: get_forgotten(vid)) {
                bag.erase(uu);
            }
        }
        assert(bag.empty());

        // frontier size
        max_frontier_size_ = 0;
        for (const auto& f: frontiers_) {
            max_frontier_size_ = std::max(static_cast<int>(f.size()), max_frontier_size_);
        }

        // map vertex to position
        // fill unused entry with V
        mapping_vertex_pos_.assign(V, V);
        std::vector<bool> used(max_frontier_size_, false);
        for (InternalVertexId vid = 0; vid < static_cast<InternalVertexId>(V); ++vid) {
            // introduce vertex
            Vertex v = get_vertex(vid);
            bool success = false;
            for (FrontierPosition i = 0; i < max_frontier_size_; ++i) {
                if (!used[i]) {
                    mapping_vertex_pos_[v] = i;
                    used[i] = true;
                    success = true;
                    break;
                }
            }
            assert(success);

            // forget vertex
            for (auto uu: get_forgotten(vid)) {
                FrontierPosition released = mapping_vertex_pos_[uu];
                used[released] = false;
            }
        }
        for (Vertex u = 0; u < static_cast<Vertex>(V); ++u) {
            assert(mapping_vertex_pos_[u] != static_cast<FrontierPosition>(V));
        }
        for (FrontierPosition i = 0; i < max_frontier_size_; ++i) {
            assert(!used[i]);
        }
    }
};

} // namespace graph
} // namespace pyzdd

#endif // PYZDD_GRAPH_H
