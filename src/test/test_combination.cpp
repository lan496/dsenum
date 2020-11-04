#include <iostream>
#include <unordered_set>
#include <cassert>

#include <tdzdd/DdSpec.hpp>
#include <tdzdd/DdStructure.hpp>

#include <type.hpp>
#include <spec/combination.hpp>

using namespace pyzdd;
using namespace pyzdd::combination;

// https://stackoverflow.com/questions/29855908/c-unordered-set-of-vectors
template<typename T>
struct VectorHash {
    size_t operator()(const std::vector<T>& v) const {
        std::hash<T> hasher;
        size_t seed = 0;
        for (auto i : v) {
            seed ^= hasher(i) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        }
        return seed;
    }
};

std::string check_enumerated(int n, int k) {
    Combination spec(n, k);
    tdzdd::DdStructure<2> dd(spec);
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

    auto expect = brute_force_combination(n, k);
    if (actual != std::to_string(expect.size())) {
        std::cerr << "The cardinality is wrong: (actual, expect) = ("
                  << actual << ", " << expect.size() << ")" << std::endl;
        exit(1);
    }
    std::unordered_set<std::vector<int>, VectorHash<int>> uset_expect;
    for (auto coloring: expect) {
        uset_expect.insert(coloring);
    }

    std::unordered_set<std::vector<int>, VectorHash<int>> uset_actual;
    for (auto itr = dd.begin(), end = dd.end(); itr != end; ++itr) {
        std::vector<int> choice;
        for (auto level: *itr) {
            choice.emplace_back(n - level);
        }
        uset_actual.insert(choice);
    }

    if (uset_actual != uset_expect) {
        std::cerr << "DD" << std::endl;
        for (auto choice: uset_actual) {
            std::cerr << "  ";
            for (auto c: choice) {
                std::cerr << c;
            }
            std::cerr << std::endl;
        }

        std::cerr << "brute force" << std::endl;
        for (auto choice: uset_expect) {
            std::cerr << "  ";
            for (auto c: choice) {
                std::cerr << c;
            }
            std::cerr << std::endl;
        }
        assert(false);
    }

    return actual;
}

void test1() {
    int n = 4;
    int k = 2;

    std::string cardinality_expect = "6";
    auto actual = check_enumerated(n, k);
    assert(actual == cardinality_expect);
}

void test_small(int n_max) {
    for (int n = 1; n <= n_max; ++n) {
        for (int k = 0; k <= n; ++k) {
            check_enumerated(n, k);
        }
    }
}

int main() {
    tdzdd::MessageHandler::showMessages(true);
    test1();
    tdzdd::MessageHandler::showMessages(false);
    test_small(8);
    return 0;
}
