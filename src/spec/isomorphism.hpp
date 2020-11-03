#ifndef PYZDD_ISOMORPHISM_H
#define PYZDD_ISOMORPHISM_H

#include <climits>
#include <cassert>
#include <type_traits>
#include <tdzdd/DdSpec.hpp>
#include "../type.hpp"
#include "../permutation.hpp"

namespace pyzdd {
namespace permutation {
namespace isomorphism {

// ============================================================================
// Internal state for DD
// ============================================================================

/// @brief manage the comparison of two colorings at a position.
enum struct ElementComp : char {
    UNKNOWN = 0,
    GREATER = 1,  // sigma(x[i]) > x[i]
    EQUAL = 2,    // sigma(x[i]) == x[i]
    SMALLER = 3,  // sigma(x[i]) < x[i]
};

enum struct CompResult : char {
    UNKNOWN = 0,
    GREATER = 1, // sigma(x) > x
    SMALLER_OR_EQUAL = 2, // sigma(x) <= x
};

// this class is POD type
struct CompPivot {
    /// the whole comparison
    CompResult result;
    /// the smallest element whose comparison is obtained
    Element pivot;
    /// comparison of (pivot, sigma(pivot))
    ElementComp comp;
};

using BinaryColor = char;  // 0 or 1
const BinaryColor UNUSED_COLOR = -1;

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

std::ostream& operator<<(std::ostream& os, const CompResult result) {
    if (result == CompResult::UNKNOWN) {
        os << " U";
    } else if (result == CompResult::GREATER) {
        os << " >";
    } else if (result == CompResult::SMALLER_OR_EQUAL){
        os << "<=";
    }
    return os;
}

std::ostream& operator<<(std::ostream& os, const CompPivot cp) {
    os << "(result='" << cp.result << "', pivot=" << static_cast<int>(cp.pivot) << ", comp='" << cp.comp << "')";
    return os;
}

std::ostream& operator<<(std::ostream& os, const BinaryColor color) {
    if (color == UNUSED_COLOR) {
        os << "U";
    } else {
        os << static_cast<int>(color);
    }
    return os;
}

// ============================================================================
// Isomorphism Elimination DD
// ============================================================================

/// @brief DD specification for representing coloring that is lexicographically
///        smaller than or equal to a permutated coloring. Note that the
///        reference enumerate lexicographically greater than or equal.
/// @details see T. Horiyama, M. Miyasaka, and R. Sasaki, in the Canadian Conference on Computational Geometry (2018).
class IsomorphismElimination:
    public tdzdd::HybridDdSpec<IsomorphismElimination, CompPivot, BinaryColor, 2> {
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
        assert(std::is_pod<CompPivot>::value);
        assert(std::is_pod<BinaryColor>::value);

        setArraySize(max_frontier_size_);
    }

    int getRoot(CompPivot& cp, BinaryColor* state) const {
        reset_comppivot(cp);
        init_state(state);
        return n_;
    }

    int getChild(CompPivot& cp, BinaryColor* state, int level, int value) const {
        const Element e = n_ - level;

        // if the whole comparison is already determined, do not care value
        if (cp.result == CompResult::SMALLER_OR_EQUAL) {
            // if level == 1, level - 1 == Terminal::REJECT
            if (level == 1) {
                return Terminal::ACCEPT;
            } else {
                return level - 1;
            }
        }
        assert(cp.result == CompResult::UNKNOWN);

        // initialize states for introduced elements
        reset_color(state, e);

#ifdef _DEBUG
        std::cerr << std::endl;
        std::cerr << "# call element=" << e << ", value=" << value << std::endl;
        std::cerr << "Before determing color" << std::endl;
        dump_comppivot_and_state(std::cerr, cp, state, level);
#endif

        // color e with value
        set_color(state, e, static_cast<BinaryColor>(value));

        // comparison
        const auto& compared = pfm_.get_compared(e);
        std::vector<std::pair<Element, ElementComp>> new_comps;
        new_comps.reserve(compared.size());
        for (auto i_and_si: compared) {
            BinaryColor color_i = get_color(state, i_and_si.first);
            BinaryColor color_si = get_color(state, i_and_si.second);
            assert(color_i != UNUSED_COLOR);
            assert(color_si != UNUSED_COLOR);

            if (color_i > color_si) {
                new_comps.emplace_back(i_and_si.first, ElementComp::GREATER);
            } else if (color_i < color_si) {
                new_comps.emplace_back(i_and_si.first, ElementComp::SMALLER);
            } else {
                new_comps.emplace_back(i_and_si.first, ElementComp::EQUAL);
            }
        }

#ifdef _DEBUG
        std::cerr << "After determing color" << std::endl;
        dump_comppivot_and_state(std::cerr, cp, state, level);
#endif

        // compress ElementComp in frontier
        compress_state(cp, state, e, new_comps);
        if (cp.result == CompResult::GREATER) {
            return Terminal::REJECT;
        }

#ifdef _DEBUG
        std::cerr << "After compression" << std::endl;
        dump_comppivot_and_state(std::cerr, cp, state, level);
#endif

        // forget
        const std::vector<Element>& forgotten = pfm_.get_forgotten(e);
        for (auto ef: forgotten) {
            // color of ef is no more needed.
            reset_color(state, ef);
        }

#ifdef _DEBUG
        std::cerr << "After forgetting element" << std::endl;
        dump_comppivot_and_state(std::cerr, cp, state, level);
#endif

        if (level == 1) {
            if (cp.result == CompResult::SMALLER_OR_EQUAL) {
                return Terminal::ACCEPT;
            } else {
                return Terminal::REJECT;
            }
        }
        return level - 1;
    }
private:
    void reset_comppivot(CompPivot& cp) const {
        cp.result = CompResult::UNKNOWN;
        cp.pivot = n_;
        cp.comp = ElementComp::EQUAL;  // set EQUAL to simplify implementation of compress_state
    }

    void init_state(BinaryColor* state) const {
        for (int i = 0; i < max_frontier_size_; ++i) {
            state[i] = UNUSED_COLOR;
        }
    }

    BinaryColor get_color(BinaryColor* state, Element e) const {
        return state[pfm_.map_element_to_position(e)];
    }

    void set_color(BinaryColor* state, Element e, BinaryColor color) const {
        state[pfm_.map_element_to_position(e)] = color;
    }

    void reset_color(BinaryColor* state, Element e) const {
        state[pfm_.map_element_to_position(e)] = UNUSED_COLOR;
    }

    /// @brief Algorithm 3 in the reference
    void compress_state(CompPivot& cp, BinaryColor* state, Element e,
                              const std::vector<std::pair<Element, ElementComp>>& new_comps) const {
        // update pivot and comp
        for (auto i_and_comp: new_comps) {
            Element ec = i_and_comp.first;
            ElementComp new_comp = i_and_comp.second;
            assert(ec != cp.pivot);
            assert(new_comp != ElementComp::UNKNOWN);
            if ((cp.pivot == static_cast<Element>(n_)) || (ec < cp.pivot)) {
                if (new_comp != ElementComp::EQUAL) {
                    // here old cp.pivot does not need to remember
                    cp.pivot = ec;
                    cp.comp = new_comp;
                }
            } else {
                if ((cp.comp == ElementComp::EQUAL) && (new_comp != ElementComp::EQUAL)) {
                    // here old cp.pivot is useless
                    cp.pivot = ec;
                    cp.comp = new_comp;
                }
            }
        }

        // update result
        const auto& comp_finished = pfm_.get_comp_finished(e);
        // if all elements less than cp.pivot is contained in comp_finished,
        // we can compare the whole colorings.
        // Since comp_finished is sorted in the ascending order,
        if ((comp_finished.size() > e) && (comp_finished[e] == e)) {
            if (comp_finished.size() == static_cast<size_t>(n_)) {
                // Here, all comparisons are finished.
                if ((cp.comp == ElementComp::SMALLER) || (cp.comp == ElementComp::EQUAL)) {
                    cp.result = CompResult::SMALLER_OR_EQUAL;
                } else if (cp.comp == ElementComp::GREATER) {
                    cp.result = CompResult::GREATER;
                } else {
                    assert(false); // unreachable
                }
            } else {
                if (cp.comp == ElementComp::SMALLER) {
                    // determined to be lexicographically smaller than the permutated.
                    reset_comppivot(cp);
                    cp.result = CompResult::SMALLER_OR_EQUAL;
                    // forget state for compression
                    init_state(state);
                } else if (cp.comp == ElementComp::GREATER) {
                    cp.result = CompResult::GREATER;
                }
            }
        }
    }

    void dump_comppivot_and_state(std::ostream& os, CompPivot& cp, BinaryColor* state, int level) const {
        Element e = n_ - level;
        os << "     CompPivot : " << cp << std::endl;

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
    }
};

bool is_lexicographically_smaller_or_equal(const std::vector<BinaryColor>& lhs, const std::vector<BinaryColor>& rhs) {
    size_t n = lhs.size();
    assert(rhs.size() == n);
    for (size_t i = 0; i < n; ++i) {
        if (lhs[i] < rhs[i]) {
            return true;
        } else if (lhs[i] > rhs[i]) {
            return false;
        }
    }
    // here, lhs == rhs
    return true;
}

/// @brief let n = perm.get_size(), this function takes O(2^n).
std::vector<std::vector<BinaryColor>> brute_force_isomophism_elimination(const Permutation& perm) {
    size_t n = perm.get_size();
    if (n > 64) {
        std::cerr << "The current implementation does not support n > 64." << std::endl;
    }
    std::vector<std::vector<BinaryColor>> winners;
    for (uint64_t bits = 0; bits < (static_cast<uint64_t>(1) << n); ++bits) {
        std::vector<BinaryColor> colors(n);
        for (size_t i = 0; i < n; ++i) {
            colors[i] = static_cast<BinaryColor>((bits >> i) & (static_cast<uint64_t>(1)));
        }
        if (is_lexicographically_smaller_or_equal(colors, perm.act(colors))) {
            winners.emplace_back(colors);
        }
    }
    return winners;
}

} // namespace isomorphism
} // namespace permutation
} // namespace pyzdd

#endif // PYZDD_ISOMORPHISM_H
