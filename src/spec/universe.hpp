#ifndef PYZDD_UNIVERSE_H
#define PYZDD_UNIVERSE_H

#include <tdzdd/DdSpec.hpp>
#include <type.hpp>

namespace pyzdd {
namespace universe {

class Universe: public tdzdd::StatelessDdSpec<Universe, 2> {
    /// number of variables
    const int n;

public:
    Universe() = delete;
    Universe(const Universe&) = default;
    Universe(const int n) : n(n) {}

    int getRoot() const {
        return n;
    }

    int getChild(int level, int value) const {
        --level;
        return (level == 0) ? Terminal::ACCEPT : level;
    }
};

} // namespace universe
} // namespace pyzdd

#endif // PYZDD_UNIVERSE_H
