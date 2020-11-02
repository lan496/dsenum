
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

} // namespace combination
} // namespace pyzdd
#endif // PYZDD_COMBINATION_H
