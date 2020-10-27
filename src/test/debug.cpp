#include <iostream>
#include <string>
#include <tdzdd/DdSpec.hpp>
#include <tdzdd/DdStructure.hpp>
#include "graph.hpp"
#include "type.hpp"
#include "spec/graphspec.hpp"
using namespace graph;
using namespace tdzdd;

Graph make_grid(int n) {
    int V = n * n;
    Graph g(V);
    for (int y = 0; y < n; ++y) {
        for (int x = 0; x < n; ++x) {
            if (y < n - 1) {
                add_undirected_edge(g, y * n + x, (y + 1) * n + x, 1);
            }
            if (x < n - 1) {
                add_undirected_edge(g, y * n + x, y * n + x + 1, 1);
            }
        }
    }
    return g;
}

Graph make_path(int n) {
    Graph g(n);
    for (int i = 0; i < n - 1; ++i) {
        add_undirected_edge(g, i, i + 1, 1);
    }
    return g;
}

Graph make_cycle(int n) {
    Graph g(n);
    for (int i = 0; i < n; ++i) {
        add_undirected_edge(g, i, (i + 1) % n, 1);
    }
    return g;
}

int main() {
    /*
    int n = 10;
    int s = 0;
    int t = n - 1;
    auto g = make_cycle(n);
    auto gaux = GraphAuxiliary(g);
    auto spec = SimpleSTPath(g, s, t);
    DdStructure<2> dd(spec);
    std::cerr << "n = " << n << ", # of solutions = " << dd.zddCardinality() << std::endl;
    for (auto itr = dd.begin(); itr != dd.end(); ++itr) {
        for (auto item: *itr) {
            std::cerr << " " << item;
        }
        std::cerr << std::endl;
    }
    */
    // https://oeis.org/A007764
    std::vector<std::string> ans = {
        "2",
        "12",
        "184",
        "8512",
        "1262816",
        "575780564",
        "789360053252",
        "3266598486981642",
        "41044208702632496804",
        "1568758030464750013214100",
        "182413291514248049241470885236",
        "64528039343270018963357185158482118",
        "69450664761521361664274701548907358996488",
    };

    for (int n = 2; n <= 10; ++n) {
        Graph g = make_grid(n);
        int s = 0;
        int t = n * n - 1;
        auto gaux = GraphAuxiliary(g);

        auto spec = SimpleSTPath(g, s, t);

        // monitor time and memory
        tdzdd::MessageHandler::showMessages(true);
        tdzdd::MessageHandler mh;
        mh.begin("begin");

        DdStructure<2> dd(spec);
        dd.zddReduce();

        mh.end();

        auto actual = dd.zddCardinality();
        auto expect = ans[n - 2];
        std::cerr << "n = " << n << ", # of solutions = " << dd.zddCardinality() << std::endl;
        // dd.dumpDot(std::cout);
        assert(actual == expect);
    }
}
