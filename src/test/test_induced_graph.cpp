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
using namespace pyzdd::graph::induced_subgrah;

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

    VertexGraphFrontierManager vgfm(g);
#ifdef _DEBUG
    vgfm.dump(std::cerr);
#endif

    Weight target = 0;
    VertexInducedSubgraphSpec spec(vgfm, target);

    tdzdd::DdStructure<2> dd(spec);

    auto expect = "6";
    auto actual = dd.zddCardinality();
#ifdef _DEBUG
    std::cerr << "# of solutions: " << actual << std::endl;
#endif
    assert(actual == expect);

#ifdef _DEBUG
    std::ofstream ofs("debug.dot");
    dd.dumpDot(ofs);
#endif

    std::vector<std::vector<bool>> enumerated_expect = {
        {0, 0, 0},
        {0, 1, 0},
        {0, 0, 1},
        {1, 0, 0},
        {1, 1, 0},
        {1, 1, 1}
    };
    std::unordered_set<std::vector<bool>> uset_expect;
    for (auto choice: enumerated_expect) {
        uset_expect.insert(choice);
    }

    std::unordered_set<std::vector<bool>> uset_actual;
    for (auto itr = dd.begin(), end = dd.end(); itr != end; ++itr) {
        uset_actual.insert(vgfm.retrieve_vertices(*itr));
    }

    assert(uset_actual == uset_expect);
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

    VertexGraphFrontierManager vgfm(g);
    // vgfm.dump(std::cerr);
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
    // test2();
    // test3();
    return 0;
}
