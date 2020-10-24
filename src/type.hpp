#ifndef PYZDD_TYPE_H
#define PYZDD_TYPE_H

namespace tdzdd {
// for each Level `l`, the (n - l)th Variable is determined
using Variable = int;
// Let `n` be the number of variables, level of the root node is `n`
using Level = int;

enum Terminal {
    ACCEPT = -1,
    REJECT = 0,
};

} // namespace tdzdd

#endif // PYZDD_TYPE_H
