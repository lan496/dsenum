#ifndef PYZDD_GRAPH_H
#define PYZDD_GRAPH_H

#include <vector>
#include <queue>
#include <set>
#include <algorithm>
#include <cassert>

namespace graph {
/*
Edge and Graph
*/
using Vertex = size_t;
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

// Internal types for variables in DD
using InternalEdgeId = int;

class GraphAuxiliary {
    // number of vertices
    size_t V;
    // number of edges
    size_t E;
    // simple undirected graph
    Graph graph;
    // Edge edge_order[i] is evaled at the i-th
    // edge_order: InternalEdgeId -> Edge
    std::vector<Edge> edge_order;
    // bags[i] is a set of vertices needed to remember just before evaling the i-th variable.
    std::vector<std::set<Vertex>> bags;
    size_t width;
public:
    GraphAuxiliary() {}
    GraphAuxiliary(const Graph& graph) : V(graph.size()), graph(graph) {
        assert(is_simple_graph(graph));

        // number of edges
        E = 0;
        for (Vertex u = 0; u < V; ++u) {
            for (auto& e: graph[u]) {
                if (e.dst > u) {
                    ++E;
                }
            }
        }

        // order edges by BFS
        edge_order.resize(E + 1);
        InternalEdgeId edge_count = 0;
        std::vector<bool> visited_vertex(V, false);
        bags.reserve(E);
        std::set<Vertex> bag;
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
                bag.insert(v);
                for (auto e: graph[v]) {
                    if (!visited_vertex[e.dst]) {
                        edge_order[edge_count++] = e;
                        // the order of `bags.emplace_back` and `bags.insert` is crucial!
                        bags.emplace_back(std::set<Vertex>(bag));
                        que.push(e.dst);
                        bag.insert(e.dst);
                    }
                }
                bag.erase(v);
            }
        }
        // sentinel
        bags.push_back(std::set<Vertex>(bag));

        assert(edge_count == E);
        assert(bag.empty());

        // width of path decomposition
        width = 0;
        for (auto& bag: bags) {
            width = std::max(width, bag.size());
        }
    }

    bool is_simple_graph(const Graph& graph) const {
        size_t V = graph.size();
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
        return V;
    }

    size_t number_of_edges() const {
        return E;
    }

    size_t get_width() const {
        return width;
    }

    const std::vector<Edge>& get_edge_order() const {
        return edge_order;
    }

    const std::vector<std::set<Vertex>>& get_bags() const {
        return bags;
    }
};

} // namespace graph

#endif // PYZDD_GRAPH_H
