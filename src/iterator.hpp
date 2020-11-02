#ifndef PYZDD_ITERATOR_H
#define PYZDD_ITERATOR_H

#include <vector>
#include <tdzdd/DdStructure.hpp>

namespace pyzdd {
    /// @brief convert selected levels in DD into selected variables.
    /// @return std::vector<bool> if the i-th value of a returned vector is true (false),
    ///         the i-th variable is selected (not selected) in the corresponding 1-path.
    std::vector<bool> variable_choice(tdzdd::DdStructure<2>::const_iterator const &, int);
} // namespace pyzdd

#endif // PYZDD_ITERATOR_H
