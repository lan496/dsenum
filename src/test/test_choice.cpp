#include <iostream>
#include <unordered_set>
#include <cassert>

#include <tdzdd/DdSpec.hpp>
#include <tdzdd/DdStructure.hpp>

#include <type.hpp>
#include <spec/choice.hpp>

using namespace pyzdd;
using namespace pyzdd::choice;

using Variables = std::vector<Variable>;

// https://stackoverflow.com/questions/29855908/c-unordered-set-of-vectors
struct VectorHash {
    size_t operator()(const std::vector<Variable>& v) const {
        std::hash<Variable> hasher;
        size_t seed = 0;
        for (auto i : v) {
            seed ^= hasher(i) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        }
        return seed;
    }
};

std::string check_enumerated(int n, int k, const std::vector<Variable>& group, bool allow_more_than) {
    Choice spec(n, k, group, allow_more_than);
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

    auto expect = brute_force_choice(n, k, group, allow_more_than);
    if (actual != std::to_string(expect.size())) {
        std::cerr << "The cardinality is wrong: (actual, expect) = ("
                  << actual << ", " << expect.size() << ")" << std::endl;
        exit(1);
    }
    std::unordered_set<Variables, VectorHash> uset_expect;
    for (auto coloring: expect) {
        uset_expect.insert(coloring);
    }

    std::unordered_set<Variables, VectorHash> uset_actual;
    for (auto itr = dd.begin(), end = dd.end(); itr != end; ++itr) {
        Variables choice;
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
    int n = 4;
    int k = 2;
    std::vector<Variable> group = {0, 1, 3};

    std::string cardinality_expect_f = "6"; // [0, 1], [0, 3], [1, 3], [0, 1, 2], [0, 2, 3], [1, 2, 3]
    auto actual_f = check_enumerated(n, k, group, false);
    assert(actual_f == cardinality_expect_f);

    std::string cardinality_expect_t = "8"; // + [0, 1, 3], [0, 1, 2, 3]
    auto actual_t = check_enumerated(n, k, group, true);
    assert(actual_t == cardinality_expect_t);
}

void test_small(int n_max) {
    for (int n = 1; n <= n_max; ++n) {
        for (int k = 0; k <= n; ++k) {
            for (int flag = 0; flag <= 1; ++flag) {
                for (uint64_t bits = 0; bits < (static_cast<uint64_t>(1) << n); ++bits) {
                    std::vector<int> group;
                    for (size_t i = 0; i < static_cast<size_t>(n); ++i) {
                        bool take = static_cast<bool>((bits >> i) & (static_cast<uint64_t>(1)));
                        if (take) {
                            group.emplace_back(i);
                        }
                    }
                    check_enumerated(n, k, group, static_cast<bool>(flag));
                }
            }
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
