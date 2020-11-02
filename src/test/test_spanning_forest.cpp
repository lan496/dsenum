#include <iostream>
#include <string>

#include <tdzdd/DdSpec.hpp>
#include <tdzdd/DdStructure.hpp>

#include "graph.hpp"
#include "type.hpp"
#include "spec/spanning_forest.hpp"

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

int main() {
    // Number of spanning trees in nxn grid graph
    // https://oeis.org/A007341
    std::vector<std::string> ans = {
        "4",
        "192",
        "100352",
        "557568000",
        "32565539635200",
        "19872369301840986112",
        "126231322912498539682594816",
        "8326627661691818545121844900397056",
        "5694319004079097795957215725765328371712000",
        "40325021721404118513276859513497679249183623593590784",
        "2954540993952788006228764987084443226815814190099484786032640000"  // n=12
    };

    tdzdd::MessageHandler::showMessages(true);
    for (int n = 2; n <= 6; ++n) {
        Graph g = make_grid(n);
        auto gaux = GraphAuxiliary(g);
        // gaux.dump(std::cerr);
        std::cerr << "Frontier size: " << gaux.get_max_frontier_size() << std::endl;

        auto spec = spanning_forest::SpanningForestSpec(gaux);

        // monitor time and memory
        tdzdd::MessageHandler mh;
        mh.begin("begin");

        DdStructure<2> dd(spec);
        dd.zddReduce();

        mh.end();

        auto actual = dd.zddCardinality();
        auto expect = ans[n - 2];
        std::cerr << "n = " << n << ", # of solutions = " << dd.zddCardinality() << std::endl;

#ifdef _DEBUG
        std::ofstream output("debug.dot");
        dd.dumpDot(output);
#endif

        assert(actual == expect);
    }
}
