#include "iterator.hpp"
#include <set>
#include "type.hpp"

namespace tdzdd {
    std::vector<bool> variable_choice(tdzdd::DdStructure<2>::const_iterator const &itr, int n) {
        // set of levels (1..=n)
        std::set<Level> selected = *itr;
        std::vector<bool> choice(n, false);
        for (auto level = selected.rbegin(); level != selected.rend(); ++level) {
            choice[n - *level] = true;
        }
        return choice;
    }

} // namespace tdzdd
