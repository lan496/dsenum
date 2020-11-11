#include <iostream>
#include <vector>
#include <string>
#include <unordered_set>
#include <cassert>

#include <tdzdd/DdSpec.hpp>
#include <tdzdd/DdStructure.hpp>

#include "type.hpp"
#include "permutation.hpp"
#include "graph.hpp"
#include "structure_enumeration.hpp"
#include "utils.hpp"

using namespace pyzdd;
using namespace pyzdd::permutation;
using namespace pyzdd::graph;
using namespace pyzdd::derivative_structure;

void check(
    int num_sites,
    int num_types,
    const tdzdd::DdStructure<2>& dd,
    const graph::VertexGraphFrontierManager& vgfm,
    const std::string cardinality_expect,
    const std::vector<std::vector<int>>& enumerated_expect)
{
    auto cardinality_actual = dd.zddCardinality();
    std::cerr << "# of structures: " << cardinality_actual << std::endl;
    if (cardinality_actual != cardinality_expect) {
        std::cerr << "The cardinality is wrong: (actual, expect) = ("
                  << cardinality_actual << ", " << cardinality_expect
                  << ")" << std::endl;
        exit(1);
    }

    std::unordered_set<std::vector<int>, VectorHash<int>> uset_expect;
    for (auto labeling: enumerated_expect) {
        uset_expect.insert(labeling);
    }
    std::unordered_set<std::vector<int>, VectorHash<int>> uset_actual;
    for (auto itr = dd.begin(), end = dd.end(); itr != end; ++itr) {
        auto labeling = vgfm.retrieve_vertices(*itr);
        uset_actual.insert(labeling);
    }
    assert(uset_actual == uset_expect);
}

void test_binary() {
    // reproduce Fig.2
    int num_sites = 4;
    int num_types = 2;
    auto c4 = Permutation(std::vector<Element>{1, 2, 3, 0});
    auto m = Permutation(std::vector<Element>{3, 2, 1, 0});
    auto automorphism = generate_group(std::vector<Permutation>{c4, m});
    assert(automorphism.size() == 8);
    auto translations = generate_group(std::vector<Permutation>{c4});
    assert(translations.size() == 4);

    int num_variables = num_sites;

    // composition constraints
    {
        tdzdd::DdStructure<2> dd;

        std::vector<std::pair<std::vector<int>, int>> composition_constraints = {
            std::make_pair(std::vector<int>{0, 1, 2, 3}, 2)
        };

        Graph cluster_graph(num_variables);
        add_undirected_edge(cluster_graph, 0, 1, 2);
        add_undirected_edge(cluster_graph, 1, 2, 2);
        add_undirected_edge(cluster_graph, 2, 3, 2);
        add_undirected_edge(cluster_graph, 3, 0, 2);
        Weight target = 2;
        VertexGraphFrontierManager vgfm(cluster_graph);

        construct_binary_derivative_structures_with_sro(
            dd,
            num_sites,
            num_types,
            automorphism,
            translations,
            composition_constraints,
            vgfm,
            target
        );

        std::string cardinality_expect = "1";
        std::vector<std::vector<int>> enumerated_expect = {
            {0, 0, 1, 1},
        };
        check(num_sites, num_types, dd, vgfm, cardinality_expect, enumerated_expect);
    }
}

int main() {
    tdzdd::MessageHandler::showMessages(true);
    test_binary();
    return 0;
}
