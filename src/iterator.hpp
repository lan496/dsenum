#ifndef PYZDD_ITERATOR_H
#define PYZDD_ITERATOR_H

#include <vector>
#include <tdzdd/DdStructure.hpp>

namespace tdzdd {
    std::vector<bool> variable_choice(tdzdd::DdStructure<2>::const_iterator const &, int);
} // namespace tdzdd

#endif // PYZDD_ITERATOR_H
