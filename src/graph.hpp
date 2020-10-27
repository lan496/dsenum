#ifndef PYZDD_GRAPH_H
#define PYZDD_GRAPH_H

#include <iostream>
#include <vector>
#include <queue>
#include <set>
#include <algorithm>
#include <cassert>
#include "type.hpp"

namespace graph {
/*
Edge and Graph
*/
using Vertex = unsigned short;
using Weight = int;
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

using Edges = std::vector<Edge>;
using Graph = std::vector<Edges>;

void add_undirected_edge(Graph &graph, Vertex u, Vertex v, Weight w) {
    graph[u].emplace_back(Edge(u, v, w));
    graph[v].emplace_back(Edge(v, u, w));
}

// Internal types for variables in DD
using InternalEdgeId = int;
// position in PodArray
using FrontierPosition = unsigned short;

// Frontier manager for ZDD represeting edge-induced sugraphs
class GraphAuxiliary {
private:
    // number of vertices
    size_t V_;
    // number of edges
    size_t E_;
    // simple undirected graph
    Graph graph_;
    // Edge edge_order[i] is processed at the i-th
    // edge_order: InternalEdgeId -> Edge
    std::vector<Edge> edge_order_;

    // just before processing the i-th edge, the states of vertices of frontiers[i] are required.
    // for each i, frontiers_[i] is sorted by ascending order.
    std::vector<std::vector<Vertex>> frontiers_;
    // when processing the i-the edge, vertices introduced[i] enter.
    std::vector<std::vector<Vertex>> introduced_;
    // after processing the i-th edge, vertices forgotten[i] are no more needed.
    std::vector<std::vector<Vertex>> forgotten_;
    // remained[i] = frontiers[i] + introduced[i] - forgotten[i] = frontiers[i+1]
    // std::vector<std::vector<Vertex>> remained_;

    // mapping_vertex_pos_[v] is position of vertex v in frontiers
    std::vector<FrontierPosition> mapping_vertex_pos_;
    // mapping_pos_vertex_: InternalEdgeId -> FrontierPosition -> Vertex
    // std::vector<std::vector<Vertex>> mapping_pos_vertex_;

    // the maximum size of frontiers
    int max_frontier_size_;

    // ************************************************************************
    // internal member functions
    // ************************************************************************
    void order_edges() {
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
                for (auto e: graph_[v]) {
                    if (!visited_vertex[e.dst]) {
                        edge_order_[edge_count++] = e;
                        que.push(e.dst);
                    }
                }
            }
        }
        assert(edge_count == static_cast<InternalEdgeId>(E));
    }

    // construct `introduced_` and `forgotten_`
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

        // map position to vertex
        /*
        mapping_pos_vertex_.resize(E);
        for (InternalEdgeId eid = 0; eid < E; ++eid) {
            size_t fs = frontiers_[eid].size();
            mapping_pos_vertex_[eid].resize(fs);
            for (FrontierPosition i = 0; i < fs; ++i) {
                mapping_pos_vertex_[eid][i] = frontiers_[eid][i];
            }
        }
        */
    }
public:
    GraphAuxiliary() {}
    GraphAuxiliary(const Graph& graph) : V_(graph.size()), graph_(graph) {
        assert(is_simple_graph(graph_));
        if (V_ > SHRT_MAX) {
            std::cerr << "The number of vertices should be smaller than " << SHRT_MAX << std::endl;
            exit(1);
        }

        // number of edges
        E_ = 0;
        for (Vertex u = 0; u < static_cast<Vertex>(V_); ++u) {
            for (auto& e: graph_[u]) {
                if (e.dst > u) {
                    ++E_;
                }
            }
        }
        // order edges by BFS
        order_edges();
        // construct frontiers, intorduced, forgotten, max_frontier_size
        construct_frontiers();
    }

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

    tdzdd::Level get_vertex_introduced_level(Vertex u) const {
        InternalEdgeId E = number_of_edges();
        for (InternalEdgeId eid = 0; eid < E; ++eid) {
            for (const Vertex &ve: get_introduced(eid)) {
                if (ve == u) return E - eid;
            }
        }
        return 0;  // isolated vertex
    }

    void print() const {
        int V = number_of_vertices();
        int E = number_of_edges();
        std::cerr << "V=" << V << std::endl;
        std::cerr << "E=" << E << std::endl;

        std::cerr << "frontiers" << std::endl;
        for (InternalEdgeId eid = 0; eid < E; ++eid) {
            std::cerr << eid << ":";
            for (const auto& u: get_frontier(eid)) {
                std::cerr << " " << u;
            }
            std::cerr << std::endl;
        }
        std::cerr << "introduced" << std::endl;
        for (InternalEdgeId eid = 0; eid < E; ++eid) {
            std::cerr << eid << ":";
            for (const auto& u: get_introduced(eid)) {
                std::cerr << " " << u;
            }
            std::cerr << std::endl;
        }
        std::cerr << "forgotten" << std::endl;
        for (InternalEdgeId eid = 0; eid < E; ++eid) {
            std::cerr << eid << ":";
            for (const auto& u: get_forgotten(eid)) {
                std::cerr << " " << u;
            }
            std::cerr << std::endl;
        }
        std::cerr << "mapping" << std::endl;
        for (Vertex u = 0; u < static_cast<Vertex>(V_); ++u) {
            std::cerr << " " << map_vertex_to_position(u);
        }
        std::cerr << std::endl;
    }
};

} // namespace graph

#endif // PYZDD_GRAPH_H
