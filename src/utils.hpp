#ifndef PYZDD_UTILS
#define PYZDD_UTILS

#include <vector>
#include <functional>

namespace pyzdd {

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

} // namespace pyzdd
#endif // PYZDD_UTILS
