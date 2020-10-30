#include <iostream>
#include <unordered_set>
#include <cassert>

#include <tdzdd/DdSpec.hpp>
#include <tdzdd/DdStructure.hpp>

#include <permutation.hpp>
#include <type.hpp>
#include <spec/isomorphism.hpp>

using namespace permutation;
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

void check_enumerated(const DdStructure<2>& dd, const Permutation& perm) {
    using Coloring = std::vector<isomorphism::BinaryColor>;
    size_t n = perm.get_size();

    auto expect = isomorphism::brute_force_isomophism_elimination(perm);
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
    } else {
        std::cerr << "consistent with brute force." << std::endl;
    }
}

void test1(bool dump_dot) {
    std::vector<Element> sigma = {2, 1, 0};
    auto perm = Permutation(sigma);
    PermutationFrontierManager pfm(perm);
    // pfm.dump(std::cerr);

    isomorphism::IsomorphismElimination spec(pfm);

    DdStructure<2> dd(spec);
    dd.zddReduce();

    auto actual = dd.zddCardinality();
    std::cout << "# of solutions: " << actual << std::endl;

    if (dump_dot) {
        std::ofstream output("debug.dot");
        dd.dumpDot(output);
    }

    assert(actual == "6");
    check_enumerated(dd, perm);
}

void test2(bool dump_dot) {
    auto fname_dot = "debug2.dot";
    auto perm = Permutation(std::vector<Element>{1, 2, 0});
    std::string cardinarlity_expect = "5";

    PermutationFrontierManager pfm(perm);
    // pfm.dump(std::cerr);

    isomorphism::IsomorphismElimination spec(pfm);

    tdzdd::MessageHandler mh;
    mh.begin("begin");

    DdStructure<2> dd(spec);
    // dd.zddReduce();

    mh.end();

    auto actual = dd.zddCardinality();
    std::cout << "# of solutions: " << actual << std::endl;

    if (dump_dot) {
        std::ofstream output(fname_dot);
        dd.dumpDot(output);
    }

    assert(actual == cardinarlity_expect);
    check_enumerated(dd, perm);
}

int main() {
    tdzdd::MessageHandler::showMessages(true);
    test1(false);
    test2(true);

    return 0;
}
