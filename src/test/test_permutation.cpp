#include <iostream>
#include <unordered_set>
#include <cassert>

#include <tdzdd/DdSpec.hpp>
#include <tdzdd/DdStructure.hpp>

#include <permutation.hpp>
#include <type.hpp>
#include <graph.hpp>
#include <spec/isomorphism.hpp>
#include <spec/spanning_forest.hpp>

using namespace permutation;
using namespace graph;
using namespace tdzdd;

// https://stackoverflow.com/questions/29855908/c-unordered-set-of-vectors
struct VectorHash {
    size_t operator()(const std::vector<isomorphism::BinaryColor>& v) const {
        std::hash<isomorphism::BinaryColor> hasher;
        size_t seed = 0;
        for (auto i : v) {
            seed ^= hasher(i) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        }
        return seed;
    }
};

std::string check_enumerated(const Permutation& perm) {
    PermutationFrontierManager pfm(perm);

    #ifdef _DEBUG
    pfm.dump(std::cerr);
    #endif


    isomorphism::IsomorphismElimination spec(pfm);
    DdStructure<2> dd(spec);
    dd.zddReduce();

    auto actual = dd.zddCardinality();
    std::cerr << "# of solutions: " << actual << std::endl;

    #ifdef _DEBUG
    std::ofstream output("debug.dot");
    dd.dumpDot(output);
    #endif

    using Coloring = std::vector<isomorphism::BinaryColor>;
    size_t n = perm.get_size();

    auto expect = isomorphism::brute_force_isomophism_elimination(perm);
    assert(actual == std::to_string(expect.size()));
    std::unordered_set<Coloring, VectorHash> uset_expect;
    for (auto coloring: expect) {
        uset_expect.insert(coloring);
    }

    std::unordered_set<Coloring, VectorHash> uset_actual;
    for (auto itr = dd.begin(), end = dd.end(); itr != end; ++itr) {
        Coloring choice(n, 0);
        for (auto level: *itr) {
            choice[n - level] = 1;
        }
        uset_actual.insert(choice);
    }

    if (uset_actual != uset_expect) {
        std::cerr << "DD" << std::endl;
        for (auto choice: uset_actual) {
            std::cerr << "  ";
            for (auto c: choice) {
                std::cerr << static_cast<int>(c);
            }
            std::cerr << std::endl;
        }

        std::cerr << "brute force" << std::endl;
        for (auto choice: uset_expect) {
            std::cerr << "  ";
            for (auto c: choice) {
                std::cerr << static_cast<int>(c);
            }
            std::cerr << std::endl;
        }
        assert(false);
    }

    return actual;
}

void test1() {
    auto perm = Permutation(std::vector<Element>{2, 1, 0});
    std::string cardinarlity_expect = "6";
    auto actual = check_enumerated(perm);
    assert(actual == cardinarlity_expect);
}

void test2() {
    auto perm = Permutation(std::vector<Element>{1, 2, 0});
    std::string cardinarlity_expect = "5";
    auto actual = check_enumerated(perm);
    assert(actual == cardinarlity_expect);
}

void test3() {
    auto perm = Permutation(std::vector<Element>{0, 2, 1});
    std::string cardinarlity_expect = "6";
    auto actual = check_enumerated(perm);
    assert(actual == cardinarlity_expect);
}

void test4() {
    auto perm = Permutation(std::vector<Element>{2, 3, 1, 0});
    std::string cardinarlity_expect = "9";
    auto actual = check_enumerated(perm);
    assert(actual == cardinarlity_expect);
}

void test_small() {
    for (int n = 1; n <= 6; ++n) {
        std::vector<Element> sigma(n);
        for (int i = 0; i < n; ++i) {
            sigma[i] = i;
        }
        do {
            auto perm = Permutation(sigma);
            perm.dump(std::cerr);
            check_enumerated(perm);
        } while (std::next_permutation(sigma.begin(), sigma.end()));
    }
}

std::vector<Permutation> generate_edge_automorphism(const std::vector<Permutation>& vertex_automorphism, const std::vector<Edge>& edge_order) {
    int E = edge_order.size();
    std::vector<Permutation> edge_automorphism;
    for (auto perm: vertex_automorphism) {
        std::vector<Element> sigma(E);
        for (int i = 0; i < E; ++i) {
            Edge e = edge_order[i];
            Edge permuted = Edge(perm.permute(e.src), perm.permute(e.dst), e.weight);
            bool success = false;
            for (int j = 0; j < E; ++j) {
                auto ej = edge_order[j];
                if (((permuted.src == ej.src) && (permuted.dst == ej.dst)) || ((permuted.src == ej.dst) && (permuted.dst == ej.src))) {
                    sigma[i] = j;
                    success = true;
                    break;
                }
            }
            assert(success);
        }
        edge_automorphism.emplace_back(Permutation(sigma));
    }
    return edge_automorphism;
}

void test_tetrahedron() {
    int V = 4;
    Graph g(V);
    add_undirected_edge(g, 0, 1, 1);
    add_undirected_edge(g, 0, 2, 1);
    add_undirected_edge(g, 0, 3, 1);
    add_undirected_edge(g, 1, 2, 1);
    add_undirected_edge(g, 1, 3, 1);
    add_undirected_edge(g, 2, 3, 1);

    // automorphism
    std::vector<Permutation> generators = {
        std::vector<Element>{0, 2, 3, 1}, // C3
        std::vector<Element>{1, 0, 3, 2}, // C2
        std::vector<Element>{0, 1, 3, 2}, // m
    };
    auto vertex_automorphism = generate_group(generators);
    assert(vertex_automorphism.size() == 24);  // isomorphic to T_h

    // automorphism on edges
    GraphAuxiliary gaux(g);
    std::vector<Edge> edge_order = gaux.get_edge_order();
    auto edge_automorphism = generate_edge_automorphism(vertex_automorphism, edge_order);

    tdzdd::MessageHandler mh;
    mh.begin("delepments of tetrahedron");

    // enumerate labeled developments
    spanning_forest::SpanningForestSpec tree_spec(gaux);
    DdStructure<2> dd(tree_spec);
    dd.zddReduce();

    auto count_labeled_actual = dd.zddCardinality();
    auto count_labeled_expect = "16";
    assert(count_labeled_actual == count_labeled_expect);

    // enumerate unlabeled developments
    std::vector<isomorphism::IsomorphismElimination> symmetry_specs;
    for (auto perm: edge_automorphism) {
        PermutationFrontierManager pfm(perm);
        isomorphism::IsomorphismElimination symmetry_spec(pfm);
        symmetry_specs.emplace_back(symmetry_spec);
    }
    for (auto spec: symmetry_specs) {
        dd.zddSubset(spec);
        dd.zddReduce();
    }

    auto count_developments_actual = dd.zddCardinality();
    auto count_developments_expect = "1";
    assert(count_developments_actual == count_developments_expect);

    mh.end();
}

void test_cube() {
    int V = 8;
    Graph g(V);
    add_undirected_edge(g, 0, 1, 1);
    add_undirected_edge(g, 1, 2, 1);
    add_undirected_edge(g, 2, 3, 1);
    add_undirected_edge(g, 3, 0, 1);
    add_undirected_edge(g, 0, 4, 1);
    add_undirected_edge(g, 1, 5, 1);
    add_undirected_edge(g, 2, 6, 1);
    add_undirected_edge(g, 3, 7, 1);
    add_undirected_edge(g, 4, 5, 1);
    add_undirected_edge(g, 5, 6, 1);
    add_undirected_edge(g, 6, 7, 1);
    add_undirected_edge(g, 7, 4, 1);

    // automorphism
    std::vector<Permutation> generators = {
        std::vector<Element>{0, 3, 7, 4, 1, 2, 6, 5}, // C3
        std::vector<Element>{1, 2, 3, 0, 5, 6, 7, 4}, // C4
        std::vector<Element>{1, 0, 3, 2, 5, 4, 7, 6}, // m
    };
    auto vertex_automorphism = generate_group(generators);
    assert(vertex_automorphism.size() == 48);  // isomorphic to O_h

    // automorphism on edges
    GraphAuxiliary gaux(g);
    std::vector<Edge> edge_order = gaux.get_edge_order();
    auto edge_automorphism = generate_edge_automorphism(vertex_automorphism, edge_order);


}

int main() {
    tdzdd::MessageHandler::showMessages(true);
    // test1();
    // test2();
    // test3();
    test4();
    // test_small();

    // test_tetrahedron();
    // test_cube();
    return 0;
}
