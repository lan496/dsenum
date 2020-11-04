
#ifndef PYZDD_COMBINATION_H
#define PYZDD_COMBINATION_H
#include <tdzdd/DdSpec.hpp>
#include "../type.hpp"

namespace pyzdd {
namespace combination {
class Combination: public tdzdd::DdSpec<Combination, int, 2> {
    const int n;
    const int k;

public:
    // choice k elements out of n.
    Combination(int n, int k) : n(n), k(k) {}

    int getRoot(int& state) const {
        state = 0;
        return n;
    }

    int getChild(int& state, int level, int value) const {
        // `state` is the number of choosen elements
        state += value;
        --level;
        if (level == 0) {
            return (state == k) ? Terminal::ACCEPT : Terminal::REJECT;
        }
        if ((state > k) || (state + level < k)) {
            return Terminal::REJECT;
        }
        return level;
    }
};

std::vector<std::vector<int>> brute_force_combination(int n, int k) {
    if (n > 64) {
        std::cerr << "The current implementation does not support n > 64." << std::endl;
    }

    std::vector<std::vector<int>> combs;
    for (uint64_t bits = 0; bits < (static_cast<uint64_t>(1) << n); ++bits) {
        std::vector<int> comb;
        int count = 0;
        for (size_t i = 0; i < static_cast<size_t>(n); ++i) {
            bool take = static_cast<bool>((bits >> i) & (static_cast<uint64_t>(1)));
            if (take) {
                comb.emplace_back(i);
                ++count;
            }
        }
        std::reverse(comb.begin(), comb.end());
        if (count == k) {
            combs.emplace_back(comb);
        }
    }
    return combs;
}

} // namespace combination
} // namespace pyzdd
#endif // PYZDD_COMBINATION_H
