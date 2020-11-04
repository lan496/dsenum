#include <iostream>
#include <vector>
#include <string>
#include <unordered_set>
#include <cassert>

#include <tdzdd/DdSpec.hpp>
// #include <tdzdd/DdStructure.hpp>

#include "type.hpp"
#include "permutation.hpp"
#include "structure_enumeration.hpp"

using namespace pyzdd;
using namespace pyzdd::permutation;
using namespace pyzdd::derivative_structure;

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

void check(
    int num_sites,
    int num_types,
    const tdzdd::DdStructure<2>& dd,
    const std::string cardinality_expect,
    const std::vector<std::vector<Element>>& enumerated_expect)
{
    auto cardinality_actual = dd.zddCardinality();
    std::cerr << "# of nonequivalent structures: " << cardinality_actual << std::endl;
    if (cardinality_actual != cardinality_expect) {
        std::cerr << "The cardinality is wrong: (actual, expect) = ("
                  << cardinality_actual << ", " << cardinality_expect
                  << ")" << std::endl;
        exit(1);
    }

    std::unordered_set<std::vector<Element>, VectorHash<Element>> uset_expect;
    for (auto labeling: enumerated_expect) {
        uset_expect.insert(labeling);
    }
    std::unordered_set<std::vector<Element>, VectorHash<Element>> uset_actual;
    for (auto itr = dd.begin(), end = dd.end(); itr != end; ++itr) {
        auto labeling = convert_to_labeling(itr, num_sites, num_types);
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

    // no option
    {
        tdzdd::DdStructure<2> dd;

        bool remove_incomplete = false;
        bool remove_superperiodic = false;
        enumerate_derivative_structures(
            num_sites,
            num_types,
            automorphism,
            translations,
            remove_incomplete,
            remove_superperiodic,
            dd
        );

        std::string cardinality_expect = "6";
        std::vector<std::vector<Element>> enumerated_expect = {
            {0, 0, 0, 0},
            {0, 0, 0, 1},
            {0, 1, 0, 1},
            {0, 0, 1, 1},
            {0, 1, 1, 1},
            {1, 1, 1, 1},
        };
        check(num_sites, num_types, dd, cardinality_expect, enumerated_expect);
    }

    // remove_incomplete
    {
        tdzdd::DdStructure<2> dd;

        bool remove_incomplete = true;
        bool remove_superperiodic = false;
        enumerate_derivative_structures(
            num_sites,
            num_types,
            automorphism,
            translations,
            remove_incomplete,
            remove_superperiodic,
            dd
        );
        std::string cardinality_expect = "4";
        std::vector<std::vector<Element>> enumerated_expect = {
            {0, 0, 0, 1},
            {0, 1, 0, 1},
            {0, 0, 1, 1},
            {0, 1, 1, 1},
        };
        check(num_sites, num_types, dd, cardinality_expect, enumerated_expect);
    }

    // remove superperiodic
    {
        tdzdd::DdStructure<2> dd;

        bool remove_incomplete = false;
        bool remove_superperiodic = true;
        enumerate_derivative_structures(
            num_sites,
            num_types,
            automorphism,
            translations,
            remove_incomplete,
            remove_superperiodic,
            dd
        );

        std::string cardinality_expect = "3";
        std::vector<std::vector<Element>> enumerated_expect = {
            {0, 0, 0, 1},
            {0, 0, 1, 1},
            {0, 1, 1, 1},
        };
        check(num_sites, num_types, dd, cardinality_expect, enumerated_expect);
    }
}

void test_multi() {
    // reproduce Fig.5(b)
    int num_sites = 4;
    int num_types = 3;
    auto c4 = Permutation(std::vector<Element>{1, 2, 3, 0});
    auto m = Permutation(std::vector<Element>{3, 2, 1, 0});
    auto automorphism = generate_group(std::vector<Permutation>{c4, m});
    assert(automorphism.size() == 8);
    auto translations = generate_group(std::vector<Permutation>{c4});
    assert(translations.size() == 4);

    // no option
    {
        tdzdd::DdStructure<2> dd;
        bool remove_incomplete = false;
        bool remove_superperiodic = false;
        enumerate_derivative_structures(
            num_sites,
            num_types,
            automorphism,
            translations,
            remove_incomplete,
            remove_superperiodic,
            dd
        );

        std::string cardinality_expect = "21";
        std::vector<std::vector<Element>> enumerated_expect = {
            {0, 0, 0, 0},
            {1, 1, 1, 1},
            {2, 2, 2, 2},
            //
            {1, 1, 1, 0},
            {1, 1, 0, 0},
            {1, 0, 1, 0},
            {1, 0, 0, 0},
            //
            {2, 2, 2, 0},
            {2, 2, 0, 0},
            {2, 0, 2, 0},
            {2, 0, 0, 0},
            //
            {2, 2, 2, 1},
            {2, 2, 1, 1},
            {2, 1, 2, 1},
            {2, 1, 1, 1},
            //
            {2, 2, 1, 0}, // (2, 2, 0, 1)
            {2, 1, 2, 0}, // (2, 0, 2, 1)
            {2, 1, 1, 0}, // (2, 0, 1, 1)
            {2, 1, 0, 1},
            {2, 0, 1, 0},
            {2, 1, 0, 0}, // (2, 0, 0, 1)
        };
        check(num_sites, num_types, dd, cardinality_expect, enumerated_expect);
    }

    // remove incomplete
    {
        tdzdd::DdStructure<2> dd;
        bool remove_incomplete = true;
        bool remove_superperiodic = false;
        enumerate_derivative_structures(
            num_sites,
            num_types,
            automorphism,
            translations,
            remove_incomplete,
            remove_superperiodic,
            dd
        );

        std::string cardinality_expect = "6";
        std::vector<std::vector<Element>> enumerated_expect = {
            {2, 2, 1, 0}, // (2, 2, 0, 1)
            {2, 1, 2, 0}, // (2, 0, 2, 1)
            {2, 1, 1, 0}, // (2, 0, 1, 1)
            {2, 1, 0, 1},
            {2, 0, 1, 0},
            {2, 1, 0, 0}, // (2, 0, 0, 1)
        };
        check(num_sites, num_types, dd, cardinality_expect, enumerated_expect);
    }

    // remove superperiodic
    {
        tdzdd::DdStructure<2> dd;
        bool remove_incomplete = false;
        bool remove_superperiodic = true;
        enumerate_derivative_structures(
            num_sites,
            num_types,
            automorphism,
            translations,
            remove_incomplete,
            remove_superperiodic,
            dd
        );

        std::string cardinality_expect = "15";
        std::vector<std::vector<Element>> enumerated_expect = {
            {1, 1, 1, 0},
            {1, 1, 0, 0},
            {1, 0, 0, 0},
            //
            {2, 2, 2, 0},
            {2, 2, 0, 0},
            {2, 0, 0, 0},
            //
            {2, 2, 2, 1},
            {2, 2, 1, 1},
            {2, 1, 1, 1},
            //
            {2, 2, 1, 0}, // (2, 2, 0, 1)
            {2, 1, 2, 0}, // (2, 0, 2, 1)
            {2, 1, 1, 0}, // (2, 0, 1, 1)
            {2, 1, 0, 1},
            {2, 0, 1, 0},
            {2, 1, 0, 0}, // (2, 0, 0, 1)
        };
        check(num_sites, num_types, dd, cardinality_expect, enumerated_expect);
    }
}

int main() {
    tdzdd::MessageHandler::showMessages(true);
    test_binary();
    test_multi();
    return 0;
}
