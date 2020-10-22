#ifndef PYZDD_UNIVERSE_H
#define PYZDD_UNIVERSE_H
#include <tdzdd/DdSpec.hpp>
#include "type.hpp"

namespace pyzdd {
class Universe : public tdzdd::StatelessDdSpec<Universe, 2> {
    const int n;

public:
    Universe(int n) : n(n) {}

    int getRoot() const {
        return n;
    }

    int getChild(int level, int value) const {
        --level;
        if (level == 0) {
            return Terminal::ACCEPT;
        }
        return level;
    }
};
} // namespace pyzdd
#endif // PYZDD_UNIVERSE_H
