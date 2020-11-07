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

using namespace pyzdd;
using namespace pyzdd::graph;
using namespace pyzdd::permutation;

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
    tdzdd::DdStructure<2> dd(spec);
    // dd.zddReduce();

    auto actual = dd.zddCardinality();
#ifdef _DEBUG
    std::cerr << "# of solutions: " << actual << std::endl;
#endif

#ifdef _DEBUG
    std::ofstream output("debug.dot");
    dd.dumpDot(output);
#endif

    using Coloring = std::vector<isomorphism::BinaryColor>;
    size_t n = perm.get_size();

    auto expect = isomorphism::brute_force_isomophism_elimination(perm);
    if (actual != std::to_string(expect.size())) {
        std::cerr << "The cardinality is wrong: (actual, expect) = ("
                  << actual << ", " << expect.size() << ")" << std::endl;
        exit(1);
    }
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

void test_small(int n_max) {
    tdzdd::MessageHandler::showMessages(false);
    for (int n = 1; n <= n_max; ++n) {
        std::vector<Element> sigma(n);
        for (int i = 0; i < n; ++i) {
            sigma[i] = i;
        }
        do {
            auto perm = Permutation(sigma);
#ifdef _DEBUG
            perm.dump(std::cerr);
#endif
            check_enumerated(perm);
        } while (std::next_permutation(sigma.begin(), sigma.end()));
    }
    tdzdd::MessageHandler::showMessages(true);
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

void test_developments(const Graph& g, const std::vector<Permutation>& vertex_automorphism,
                       std::string count_labeled_expect, std::string count_developments_expect) {
    // automorphism on edges
    GraphAuxiliary gaux(g);
    std::vector<Edge> edge_order = gaux.get_edge_order();
    auto edge_automorphism = generate_edge_automorphism(vertex_automorphism, edge_order);

    tdzdd::MessageHandler mh;
    mh.begin("delepments");

    // enumerate labeled developments
    spanning_forest::SpanningForestSpec tree_spec(gaux);
    tdzdd::DdStructure<2> dd(tree_spec);
    dd.zddReduce();

    auto count_labeled_actual = dd.zddCardinality();
    assert(count_labeled_actual == count_labeled_expect);

    // sort permutations by max frontier sizes
    std::sort(edge_automorphism.begin(), edge_automorphism.end(),
              [](const Permutation& lhs, const Permutation& rhs) {
                  auto size_lhs = PermutationFrontierManager(lhs).get_max_frontier_size();
                  auto size_rhs = PermutationFrontierManager(rhs).get_max_frontier_size();
                  return size_lhs < size_rhs;
              });

    // enumerate unlabeled developments
    std::vector<isomorphism::IsomorphismElimination> symmetry_specs;
    symmetry_specs.reserve(edge_automorphism.size());
    for (auto perm: edge_automorphism) {
        PermutationFrontierManager pfm(perm);
        isomorphism::IsomorphismElimination symmetry_spec(pfm);
        symmetry_specs.emplace_back(symmetry_spec);
    }

    for (auto spec: symmetry_specs) {
        dd.zddSubset(spec);
        dd.zddReduce();
    }
#ifdef _DEBUG
    std::ofstream output("debug.dot");
    dd.dumpDot(output);
#endif

    auto count_developments_actual = dd.zddCardinality();
    if (count_developments_actual != count_developments_expect) {
        std::cerr << "The cardinality of unlabeled developments is wrong: (actual, expect) = ("
                  << count_developments_actual << ", " << count_developments_expect << ")" << std::endl;
        exit(1);
    }
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

    auto count_labeled_expect = "16";
    auto count_developments_expect = "2";
    test_developments(g, vertex_automorphism, count_labeled_expect, count_developments_expect);
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

    size_t group_order_expect = 48;
    assert(vertex_automorphism.size() == group_order_expect);  // isomorphic to O_h

    auto count_labeled_expect = "384";
    auto count_developments_expect = "11";
    test_developments(g, vertex_automorphism, count_labeled_expect, count_developments_expect);
}

void test_dodecahedron() {
    int V = 20;
    Graph g(V);
    add_undirected_edge(g, 0, 1, 1);
    add_undirected_edge(g, 0, 4, 1);
    add_undirected_edge(g, 0, 5, 1);
    add_undirected_edge(g, 1, 2, 1);
    add_undirected_edge(g, 1, 6, 1);
    add_undirected_edge(g, 2, 3, 1);
    add_undirected_edge(g, 2, 7, 1);
    add_undirected_edge(g, 3, 4, 1);
    add_undirected_edge(g, 3, 8, 1);
    add_undirected_edge(g, 4, 9, 1);
    add_undirected_edge(g, 5, 10, 1);
    add_undirected_edge(g, 5, 14, 1);
    add_undirected_edge(g, 6, 10, 1);
    add_undirected_edge(g, 6, 11, 1);
    add_undirected_edge(g, 7, 11, 1);
    add_undirected_edge(g, 7, 12, 1);
    add_undirected_edge(g, 8, 12, 1);
    add_undirected_edge(g, 8, 13, 1);
    add_undirected_edge(g, 9, 13, 1);
    add_undirected_edge(g, 9, 14, 1);
    add_undirected_edge(g, 10, 15, 1);
    add_undirected_edge(g, 11, 16, 1);
    add_undirected_edge(g, 12, 17, 1);
    add_undirected_edge(g, 13, 18, 1);
    add_undirected_edge(g, 14, 19, 1);
    add_undirected_edge(g, 15, 16, 1);
    add_undirected_edge(g, 15, 19, 1);
    add_undirected_edge(g, 16, 17, 1);
    add_undirected_edge(g, 17, 18, 1);
    add_undirected_edge(g, 18, 19, 1);

    // automorphism
    std::vector<Permutation> generators = {
        std::vector<Element>{1, 2, 3, 4, 0,
                             6, 7, 8, 9, 5,
                             11, 12, 13, 14, 10,
                             16, 17, 18, 19, 15}, // C5
        std::vector<Element>{0, 4, 3, 2, 1,
                             5, 9, 8, 7, 6,
                             14, 13, 12, 11, 10,
                             19, 18, 17, 16, 15}, // m
        std::vector<Element>{0, 4, 9, 14, 5, 1, 3, 13, 19, 10, 2, 8, 18, 15, 6, 7, 12, 17, 16, 11}, // C3
    };
    auto vertex_automorphism = generate_group(generators);  // isomorphic to A5 x Z2

    size_t group_order_expect = 120;
    assert(vertex_automorphism.size() == group_order_expect);

    auto count_labeled_expect = "5184000";
    auto count_developments_expect = "43380";
    test_developments(g, vertex_automorphism, count_labeled_expect, count_developments_expect);
}

int main() {
    tdzdd::MessageHandler::showMessages(true);
    test1();
    test2();
    test3();
    test4();
    test_small(6);
    test_tetrahedron();
    test_cube();
    test_dodecahedron();
    return 0;
}
