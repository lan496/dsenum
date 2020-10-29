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
    GREATER_OR_EQUAL = 1,  // sigma(x) >= x
    SMALLER = 2,  // sigma(x) < x
};

std::ostream& operator<<(std::ostream& os, const CompResult result) {
    if (result == CompResult::UNKNOWN) {
        os << "U";
    } else if (result == CompResult::GREATER_OR_EQUAL) {
        os << ">=";
    } else if (result == CompResult::SMALLER){
        os << "<";
    }
    return os;
}

enum struct ElementComp : char {
    UNKNOWN = 0,
    GREATER = 1,  // sigma(x[i]) > x[i]
    EQUAL = 2,    // sigma(x[i]) == x[i]
    SMALLER = 3,  // sigma(x[i]) < x[i]
};

std::ostream& operator<<(std::ostream& os, const ElementComp comp) {
    if (comp == ElementComp::UNKNOWN) {
        os << "U";
    } else if (comp == ElementComp::GREATER) {
        os << ">";
    } else if (comp == ElementComp::EQUAL) {
        os << "=";
    } else if (comp == ElementComp::SMALLER) {
        os << "<";
    }
    return os;
}

using BinaryColor = char;  // 0 or 1
const BinaryColor UNUSED_COLOR = -1;

std::ostream& operator<<(std::ostream& os, const BinaryColor color) {
    if (color == UNUSED_COLOR) {
        os << "U";
    } else {
        os << static_cast<int>(color);
    }
    return os;
}

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

        // if result == 1, do not care value
        if (result == CompResult::GREATER_OR_EQUAL) {
            if (level == 1) {
                return tdzdd::Terminal::ACCEPT;
            } else {
                return level - 1;
            }
        }

        // initialize states for introduced elements
        const std::vector<Element>& introduced = pfm_.get_introduced(e);
        for (auto ei: introduced) {
            set_color(state, ei, UNUSED_COLOR);
            set_comp(state, ei, ElementComp::UNKNOWN);
        }

        std::cerr << std::endl;
        std::cerr << "# call element=" << e << ", value=" << value << std::endl;
        std::cerr << "before processing element" << std::endl;
        dump_result_and_state(std::cerr, result, state, level);

        // color e with value
        set_color(state, e, static_cast<BinaryColor>(value));

        // comparison
        const auto& compared = pfm_.get_compared(e);
        for (auto i_and_si: compared) {
            BinaryColor color_i = get_color(state, i_and_si.first);
            BinaryColor color_si = get_color(state, i_and_si.second);
            assert(color_i != UNUSED_COLOR);
            assert(color_si != UNUSED_COLOR);
            assert(get_comp(state, i_and_si.first) == ElementComp::UNKNOWN);

            if (color_i > color_si) {
                set_comp(state, i_and_si.first, ElementComp::GREATER);
            } else if (color_i < color_si) {
                set_comp(state, i_and_si.first, ElementComp::SMALLER);
            } else {
                set_comp(state, i_and_si.first, ElementComp::EQUAL);
            }
        }

        std::cerr << "after comparison" << std::endl;
        dump_result_and_state(std::cerr, result, state, level);

        // compress ElementComp in frontier
        CompResult result_tmp = compress_state(state, e);
        if (result_tmp == CompResult::GREATER_OR_EQUAL) {
            result = CompResult::GREATER_OR_EQUAL;
            // erase state
            const auto& frontier = pfm_.get_frontier(e);
            for (auto efr: frontier) {
                set_color(state, efr, UNUSED_COLOR);
                set_comp(state, efr, ElementComp::UNKNOWN);
            }
        } else if (result_tmp == CompResult::SMALLER) {
            return tdzdd::Terminal::REJECT;
        }
        std::cerr << "after compression" << std::endl;
        dump_result_and_state(std::cerr, result, state, level);

        // forget
        const std::vector<Element>& forgotten = pfm_.get_forgotten(e);
        for (auto ef: forgotten) {
            set_color(state, ef, UNUSED_COLOR);
            set_comp(state, ef, ElementComp::UNKNOWN);
        }

        std::cerr << "after processing element" << std::endl;
        dump_result_and_state(std::cerr, result, state, level);

        if (level == 1) {
            if (result == CompResult::GREATER_OR_EQUAL) {
                return tdzdd::Terminal::ACCEPT;
            } else {
                return tdzdd::Terminal::REJECT;
            }
        }
        return level - 1;
    }
private:
    BinaryColor get_color(FrontierData* state, Element e) const {
        return state[pfm_.map_element_to_position(e)].color;
    }

    void set_color(FrontierData* state, Element e, BinaryColor color) const {
        state[pfm_.map_element_to_position(e)].color = color;
    }

    ElementComp get_comp(FrontierData* state, Element e) const {
        return state[pfm_.map_element_to_position(e)].comp;
    }

    void set_comp(FrontierData* state, Element e, ElementComp comp) const {
        state[pfm_.map_element_to_position(e)].comp = comp;
    }

    CompResult compress_state(FrontierData* state, Element e) const {
        const auto& frontier = pfm_.get_frontier(e);  // ascending order
        bool is_same = true;
        for (Element efr: frontier) {
            auto comp_efr = get_comp(state, efr);
            if (comp_efr == ElementComp::GREATER) {
                return CompResult::GREATER_OR_EQUAL;
            } else if (comp_efr == ElementComp::SMALLER) {
                return CompResult::SMALLER;
            } else if (comp_efr == ElementComp::UNKNOWN){
                // cannot judge the whole result of comparison
                is_same = false;
                break;
            }
        }
        if (is_same) {
            // all comp are equal
            return CompResult::GREATER_OR_EQUAL;
        }

        // here ['=',..., '=', 'UNKNOWN',......]
        bool need_to_remember = true;
        for (Element efr: frontier) {
            auto comp_efr = get_comp(state, efr);
            if (need_to_remember && ((comp_efr == ElementComp::GREATER) || (comp_efr == ElementComp::SMALLER))) {
                need_to_remember = false;
            } else if (!need_to_remember) {
                set_comp(state, efr, ElementComp::UNKNOWN);
                // color is already erased
            }
        }
        return CompResult::UNKNOWN;
    }

    void dump_result_and_state(std::ostream& os, CompResult& result, FrontierData* state, int level) const {
        Element e = n_ - level;
        os << "     result    : " << result << std::endl;

        os << "     frontier  :";
        const auto& frontier = pfm_.get_frontier(e);
        for (auto efr: frontier) {
            os << " " << static_cast<int>(efr);
        }
        os << std::endl;

        os << "     color     :";
        for (auto efr: frontier) {
            os << " " << get_color(state, efr);
        }
        os << std::endl;

        os << "     comparison:";
        for (auto efr: frontier) {
            os << " " << get_comp(state, efr);
        }
        os << std::endl;
    }
};

} // namespace permutation

#endif // PYZDD_ISOMORPHISM_H
