#ifndef PYZDD_GRAPH_H
#define PYZDD_GRAPH_H

#include <vector>
#include <queue>
#include <set>
#include <cassert>

namespace graph {
/*
Edge and Graph
*/
using Vertex = int;
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
    int V;
    // number of edges
    int E;
    // simple undirected graph
    Graph graph;
    // Edge edge_order[i] is evaled at the i-th
    // edge_order: InternalEdgeId -> Edge
    std::vector<Edge> edge_order;
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
        edge_order.resize(E);
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
                for (auto e: graph[u]) {
                    edge_order[edge_count++] = e;
                    if (!visited_vertex[e.dst]) {
                        que.push(e.dst);
                    }
                }
            }
        }
        assert(edge_count == E);
    }

    bool is_simple_graph(const Graph& graph) {
        int V = graph.size();
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
}

} // namespace graph

#endif PYZDD_GRAPH_H
