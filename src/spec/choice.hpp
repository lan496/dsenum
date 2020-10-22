#ifndef PYZDD_CHOISE_H
#define PYZDD_CHOISE_H
#include <vector>
#include <algorithm>
#include <tdzdd/DdSpec.hpp>
#include <tdzdd/DdStructure.hpp>
#include "type.hpp"

namespace pyzdd {
class Choice: public tdzdd::DdSpec<Choice, int, 2> {
public:
    int n;  // number of variables
    int k;  // number of choosen variables
    std::vector<Variable> group;  // list of 0-indexed variables

    Choice(int n, int k, const std::vector<int>& group) : n(n), k(k), group(group) {
        std::sort(group.begin(), group.end());
    }

    int getRoot(int& state) const {
        state = 0;
        return n;
    }

    int getChild(int& state, int level, int value) const {
        Variable idx = n - level;
        auto bounds = std::equal_range(group.begin(), group.end(), idx);
        if (bounds.second - bounds.first != 0) {
            state += value;
        }

        --level;
        if (level == 0) return (state == k) ? Terminal::ACCEPT : Terminal::REJECT;
        if (state > k) return Terminal::REJECT;
        if (state + (group.end() - bounds.second) < k) return Terminal::REJECT;
        return level;
    }
};
} // namespace pyzdd
#endif // PYZDD_CHOISE_H
