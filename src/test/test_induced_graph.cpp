#include <iostream>
#include <string>
#include <unordered_set>

#include <tdzdd/DdSpec.hpp>
#include <tdzdd/DdStructure.hpp>
#include "graph.hpp"
#include "type.hpp"
#include "utils.hpp"
#include "spec/induced_subgraph.hpp"

using namespace pyzdd;
using namespace pyzdd::graph;
using namespace pyzdd::graph::induced_subgraph;

void check(const Graph& g, Weight target, std::string cardinality, const std::vector<std::vector<int>>& enumerated_expect, bool debug) {
    VertexGraphFrontierManager vgfm(g);
    if (debug) {
        vgfm.dump(std::cerr);
    }

    VertexInducedSubgraphSpec spec(vgfm, target);

    tdzdd::DdStructure<2> dd(spec);

    auto actual = dd.zddCardinality();
    if (debug) {
        std::cerr << "# of solutions: " << actual << std::endl;
        std::ofstream ofs("debug.dot");
        dd.dumpDot(ofs);
    }
    assert(actual == cardinality);

    std::unordered_set<std::vector<int>, VectorHash<int>> uset_expect;
    for (auto choice: enumerated_expect) {
        uset_expect.insert(choice);
    }

    std::unordered_set<std::vector<int>, VectorHash<int>> uset_actual;
    for (auto itr = dd.begin(), end = dd.end(); itr != end; ++itr) {
        uset_actual.insert(vgfm.retrieve_vertices(*itr));
    }

    assert(uset_actual == uset_expect);
}

void test1() {
    /*
          1   -1
        o---o---o
        v0  v2  v1
    */
    int V = 3;
    Graph g(V);
    add_undirected_edge(g, 0, 2, 1);
    add_undirected_edge(g, 2, 1, -1);
    Weight target = 0;


    auto cardinality = "6";
    std::vector<std::vector<int>> enumerated_expect = {
        {0, 0, 0},
        {0, 1, 0},
        {0, 0, 1},
        {1, 0, 0},
        {1, 1, 0},
        {1, 1, 1}
    };
    check(g, target, cardinality, enumerated_expect, false);
}

void test2() {
    /*
       v0    v2
        o----o
        |    |
        o----o
       v3    v1
    */
    int V = 4;
    Graph g(V);
    add_undirected_edge(g, 0, 2, 1);
    add_undirected_edge(g, 0, 3, 1);
    add_undirected_edge(g, 2, 1, 1);
    add_undirected_edge(g, 3, 1, 1);
    Weight target = 2;

    auto cardinality = "4";
    std::vector<std::vector<int>> enumerated_expect = {
        {0, 1, 1, 1},
        {1, 0, 1, 1},
        {1, 1, 0, 1},
        {1, 1, 1, 0}
    };
    check(g, target, cardinality, enumerated_expect, false);
}

void test3() {
    /*
       v0  v2  v4
        o   o   o
        |   |
        o   o
       v1   v3
    */
    int V = 5;
    Graph g(V);
    add_undirected_edge(g, 0, 1, 1);
    add_undirected_edge(g, 2, 3, 1);

    VertexGraphFrontierManager vgfm(g);
    vgfm.dump(std::cerr);
}

int main() {
    test1();
    test2();
    // test3();
    return 0;
}
