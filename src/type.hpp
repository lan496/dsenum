#ifndef PYZDD_TYPE_H
#define PYZDD_TYPE_H

namespace pyzdd {
// Let `n` be the number of variables, level of the root node is `n`
using Level = int;

/// @brief 1-terminal and 0-terminal nodes
enum Terminal {
    ACCEPT = -1,
    REJECT = 0,
};

} // namespace pyzdd

#endif // PYZDD_TYPE_H
