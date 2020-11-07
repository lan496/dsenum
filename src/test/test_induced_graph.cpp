#include <iostream>
#include <string>
#include <tdzdd/DdSpec.hpp>
#include <tdzdd/DdStructure.hpp>
#include "graph.hpp"
#include "type.hpp"

using namespace pyzdd;
using namespace pyzdd::graph;

void test1() {
    /*
          1   -2
        o---o---o
        v0  v2  v1
    */
    int V = 3;
    Graph g(V);
    add_undirected_edge(g, 0, 2, 1);
    add_undirected_edge(g, 2, 1, -2);

    VertexGraphFrontierManager vgfm(g);
    // vgfm.dump(std::cerr);
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
    test2();
    test3();
    return 0;
}
