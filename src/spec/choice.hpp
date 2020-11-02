#ifndef PYZDD_CHOISE_H
#define PYZDD_CHOISE_H

#include <vector>
#include <algorithm>
#include <limits>

#include <tdzdd/DdSpec.hpp>
#include <tdzdd/DdStructure.hpp>
#include "../type.hpp"

namespace pyzdd {
namespace choice {
// for each Level `l`, the (n - l)th Variable is determined
using Variable = int;

/// @brief DD specification for k-combinations out of n elements
class Choice: public tdzdd::DdSpec<Choice, int, 2> {
    /// number of variables
    int n;
    /// number of selected variables
    int k;
    /// list of 0-indexed variables. sorted in the ascending order.
    std::vector<Variable> group;
    /// if true, allow to take more than k
    bool allow_more_than;

public:
    Choice() = default;
    Choice(const Choice&) = default;

    Choice(int n, int k, const std::vector<int>& v, bool allow_more_than) :
        n(n),
        k(k),
        group(v),
        allow_more_than(allow_more_than)
    {
        if (n > std::numeric_limits<int>::max()) {
            std::cerr << "The number of vertices should be smaller than "
                      << std::numeric_limits<int>::max() << std::endl;
            exit(1);
        }

        for (auto vv: v) {
            if ((vv < 0) || (vv >= n)) {
                std::cerr << "v has an invalid variable: " << vv << std::endl;
                exit(1);
            }
        }

        std::sort(group.begin(), group.end());
    }

    int getRoot(int& state) const {
        state = 0;
        return n;
    }

    int getChild(int& state, Level level, int value) const {
        Variable idx = n - level;
        auto bounds = std::equal_range(group.begin(), group.end(), idx);
        if (bounds.second - bounds.first != 0) {
            // once the number of choices exceeds k, we do not need to count the excesses
            state = std::min(state + value, k + 1);
        }

        --level;
        if (level == 0) return satisfy(state) ? Terminal::ACCEPT : Terminal::REJECT;
        if (!allow_more_than && (state > k)) return Terminal::REJECT;
        if (state + (group.end() - bounds.second) < k) return Terminal::REJECT;
        return level;
    }
private:
    bool satisfy(int state) const {
        if (allow_more_than) {
            return state >= k;
        } else {
            return state == k;
        }
    }
};

std::vector<std::vector<Variable>> brute_force_choice(int n,
                                                      int k,
                                                      const std::vector<Variable>& group,
                                                      bool allow_more_than)
{
    if (n > 64) {
        std::cerr << "The current implementation does not support n > 64." << std::endl;
    }

    std::vector<std::vector<Variable>> combs;
    for (uint64_t bits = 0; bits < (static_cast<uint64_t>(1) << n); ++bits) {
        std::vector<Variable> comb;
        int count = 0;
        for (size_t i = 0; i < static_cast<size_t>(n); ++i) {
            bool take = static_cast<bool>((bits >> i) & (static_cast<uint64_t>(1)));
            if (take) {
                comb.emplace_back(i);
                auto bounds = std::equal_range(group.begin(), group.end(), i);
                if (bounds.second - bounds.first != 0) {
                    count += 1;
                }
            }
        }
        std::reverse(comb.begin(), comb.end());
        if (count == k) {
            combs.emplace_back(comb);
        } else if ((count > k) && allow_more_than) {
            combs.emplace_back(comb);
        }
    }
    return combs;
}

} // namespace choice
} // namespace pyzdd
#endif // PYZDD_CHOISE_H
