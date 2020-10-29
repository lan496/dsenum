#ifndef PYZDD_ISOMORPHISM_H
#define PYZDD_ISOMORPHISM_H

#include <climits>
#include <cassert>
#include <type_traits>
#include <tdzdd/DdSpec.hpp>
#include "../type.hpp"
#include "../permutation.hpp"

namespace permutation {

enum struct CompResult : char {
    UNKNOWN = 0,
    GREATER = 1,  // sigma(x) > x
    EQUAL = 2,    // sigma(x) == x
    SMALLER = 3,  // sigma(x) < x
};

enum struct ElementComp : char {
    UNKNOWN = 0,
    GREATER = 1,  // sigma(x[i]) > x[i]
    EQUAL = 2,    // sigma(x[i]) == x[i]
    SMALLER = 3,  // sigma(x[i]) < x[i]
};

using BinaryColor = char;  // 0 or 1
const BinaryColor UNUSED_COLOR = -1;

struct FrontierData {
    BinaryColor color;
    ElementComp comp;
};

class IsomorphismElimination:
    public tdzdd::HybridDdSpec<IsomorphismElimination, CompResult, FrontierData, 2> {
    const PermutationFrontierManager pfm_;
    const int max_frontier_size_;
    // the size of permutation
    const int n_;
public:
    IsomorphismElimination() = delete;
    IsomorphismElimination(const PermutationFrontierManager& pfm) :
        pfm_(pfm),
        max_frontier_size_(pfm.get_max_frontier_size()),
        n_(pfm.get_size())
    {
        // sanity check on types
        assert(std::is_pod<FrontierData>::value);

        setArraySize(max_frontier_size_);
    }

    int getRoot(CompResult& result, FrontierData* state) const {
        result = CompResult::UNKNOWN;
        for (int i = 0; i < max_frontier_size_; ++i) {
            state[i].color = UNUSED_COLOR;
            state[i].comp = ElementComp::UNKNOWN;
        }
        return n_;
    }

    int getChild(CompResult& result, FrontierData* state, int level, int value) const {
        const Element e = n_ - level;

        // TODO
        // initialize states for introduced elements

        // color e with value

        // compare with e

        // forget and branch on determined elements

        // if (level == 1) return tdzdd::Terminal::ACCEPT;
        return level - 1;
    }
private:
    BinaryColor get_color(FrontierData* state, Element e) const {
        return state[pfm_.map_element_to_position(e)].color;
    }

    void set_color(FrontierData* state, Element e, BinaryColor color) {
        state[pfm_.map_element_to_position(e)].color = color;
    }

    ElementComp get_comp(FrontierData* state, Element e) const {
        return state[pfm_.map_element_to_position(e)].comp;
    }

    void set_comp(FrontierData* state, Element e, ElementComp comp) {
        state[pfm_.map_element_to_position(e)].comp = comp;
    }
};

} // namespace permutation

#endif // PYZDD_ISOMORPHISM_H
