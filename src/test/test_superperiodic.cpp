#include <iostream>
#include <unordered_set>
#include <cassert>

#include <tdzdd/DdSpec.hpp>
#include <tdzdd/DdStructure.hpp>

#include <permutation.hpp>
#include <type.hpp>
#include <graph.hpp>
#include <spec/superperiodic.hpp>

using namespace permutation;
using namespace tdzdd;

// https://stackoverflow.com/questions/29855908/c-unordered-set-of-vectors
struct VectorHash {
    size_t operator()(const std::vector<superperiodic::BinaryColor>& v) const {
        std::hash<superperiodic::BinaryColor> hasher;
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

    superperiodic::SuperperiodicElimination spec(pfm);
    DdStructure<2> dd(spec);
#ifndef _DEBUG
    dd.zddReduce();
#endif

    auto actual = dd.zddCardinality();
#ifdef _DEBUG
    std::cerr << "# of solutions: " << actual << std::endl;
#endif

    #ifdef _DEBUG
    std::ofstream output("debug.dot");
    dd.dumpDot(output);
    #endif

    using Coloring = std::vector<superperiodic::BinaryColor>;
    size_t n = perm.get_size();

    auto expect = superperiodic::brute_force_superperiodic_elimination(perm);
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
    std::string cardinality_expect = "4"; // [1, 0, 0], [1, 1, 0], [0, 0, 1], [0, 1, 1]
    auto actual = check_enumerated(perm);
    assert(actual == cardinality_expect);
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

int main() {
    tdzdd::MessageHandler::showMessages(true);
    test1();
    test_small(6);
    return 0;
}
