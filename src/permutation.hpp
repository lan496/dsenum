#ifndef PYZDD_PERMUTATION_H
#define PYZDD_PERMUTATION_H

#include <iostream>
#include <algorithm>
#include <vector>
#include <set>
#include <limits>
#include <cassert>
#include <type_traits>
#include "type.hpp"

namespace permutation {
// unsigned
using Element = unsigned int;
// position of Element in PodArray
using FrontierPosition = Element;

// ============================================================================
// Permutation class
// ============================================================================
/*
# Arguments
- sigma_:
    represents a permutation in "one-line" notation.
    That is, `sigma` moves `i` to `sigma[i]`.
*/
class Permutation {
    const Element n_;
    const std::vector<Element> sigma_;
public:
    Permutation() = delete;
    Permutation(const std::vector<Element>& sigma) :
        n_(static_cast<Element>(sigma.size())),
        sigma_(sigma)
    {
        // sanity check on Element class
        assert(std::is_unsigned<Element>::value);
        if (n_ > std::numeric_limits<Element>::max()) {
            std::cerr << "The number of elements in a permutation should be less than "
                        << std::numeric_limits<Element>::max() << "." << std::endl;
        }

        if (!is_permutation(sigma)) {
            std::cerr << "Given list is not a permutation!" << std::endl;
            exit(1);
        }
    }

    size_t get_size() const {
        return n_;
    }

    Element permute(const Element i) const {
        return sigma_[i];
    }

    template<typename COLOR>
    std::vector<COLOR> act(const std::vector<COLOR> colors) const {
        std::vector<COLOR> permutated(colors.size());
        for (Element i = 0, n = get_size(); i < n; ++i) {
            permutated[permute(i)] = colors[i];
        }
        return permutated;
    }

    void dump(std::ostream& os) const {
        os << "(";
        for (Element i = 0, n = get_size(); i < n; ++i) {
            os << " " << permute(i);
        }
        os << " )";
    }

private:
    bool is_permutation(const std::vector<Element>& sigma) const {
        auto n = get_size();
        std::vector<bool> visited(false, n);
        for (size_t i = 0; i < n; ++i) {
            if ((sigma[i] > n) || visited[sigma[i]]) {
                return false;
            }
        }
        return true;
    }
};

// ============================================================================
// Prepare frontiers for Isomophism Elimination
// ============================================================================
class PermutationFrontierManager {
    // permutation
    Permutation perm_;
    // just before coloring the i-th element, the states of elements in
    // frontiers_[i] are required.
    std::vector<std::vector<Element>> frontiers_;
    // when processing the i-th element, elements in introduced_[i] enter.
    std::vector<std::vector<Element>> introduced_;
    // just after coloring the i-th element, elements in compared_[i] can be
    // compared with the i-th color.
    std::vector<std::vector<Element>> compared_;
    // after processing the i-th element, elements in forgotten_[i] are no more
    // required.
    std::vector<std::vector<Element>> forgotten_;
    // mapping_element_to_pos_[e] is a potision of element e in frontiers.
    // Be careful several elements may share the same position!
    std::vector<FrontierPosition> mapping_element_to_pos_;
    // the maximum size of frontiers
    int max_frontier_size_;
public:
    PermutationFrontierManager() = delete;
    PermutationFrontierManager(const Permutation& perm) : perm_(perm) {
        // construct frontiers, forgotten, and max_frontier_size
        construct();
    }

    Element get_size() const {
        return perm_.get_size();
    }

    int get_max_frontier_size() const {
        return max_frontier_size_;
    }

    const std::vector<Element>& get_frontier(Element e) const {
        return frontiers_[e];
    }

    const std::vector<Element>& get_introduced(Element e) const {
        return introduced_[e];
    }

    const std::vector<Element>& get_compared(Element e) const {
        return compared_[e];
    }

    const std::vector<Element>& get_forgotten(Element e) const {
        return forgotten_[e];
    }

    const FrontierPosition map_element_to_position(Element e) const {
        return mapping_element_to_pos_[e];
    }

    void dump(std::ostream& os) const {
        // show permutation
        os << "permutation: ";
        perm_.dump(os);
        os << std::endl;

        // show frontier related
        auto n = get_size();

        os << "frontiers" << std::endl;
        for (Element e = 0; e < n; ++e) {
            os << e << ":";
            for (auto ef: get_frontier(e)) {
                os << " " << ef;
            }
            os << std::endl;
        }

        os << "compared" << std::endl;
        for (Element e = 0; e < n; ++e) {
            os << e << ":";
            for (auto ef: get_compared(e)) {
                os << " " << ef;
            }
            os << std::endl;
        }

        os << "forgotton" << std::endl;
        for (Element e = 0; e < n; ++e) {
            os << e << ":";
            for (auto ef: get_forgotten(e)) {
                os << " " << ef;
            }
            os << std::endl;
        }

        os << "mapping" << std::endl;
        for (Element e = 0; e < n; ++e) {
            os << " " << map_element_to_position(e);
            os << std::endl;
        }
    }
private:
    void construct() {
        auto n = get_size();

        // prepare introduced_
        introduced_.resize(n);
        std::vector<Element> visited_element(n, false);
        for (Element e = 0; e < n; ++e) {
            if (!visited_element[e]) {
                introduced_[e].emplace_back(e);
                visited_element[e] = true;
            }
            Element e_perm = perm_.permute(e);
            if (!visited_element[e_perm]) {
                introduced_[e].emplace_back(e_perm);
                visited_element[e_perm] = true;
            }
        }
        for (Element e = 0; e < n; ++e) {
            assert(visited_element[e]);
        }

        // prepare forgotten_ by backward
        std::vector<bool> processed_element(n, false);
        for (int e = n - 1; e >= 0; --e) {  // Be careful overflow!
            if (!processed_element[static_cast<Element>(e)]) {
                forgotten_[static_cast<Element>(e)].emplace_back(static_cast<Element>(e));
                processed_element[static_cast<Element>(e)] = true;
            }
            Element e_perm = perm_.permute(static_cast<Element>(e));
            if (!processed_element[e_perm]) {
                forgotten_[static_cast<Element>(e)].emplace_back(e_perm);
                processed_element[e_perm] = true;
            }
        }
        for (Element e = 0; e < n; ++e) {
            assert(processed_element[e]);
        }

        // prepare frontiers
        frontiers_.resize(n);
        std::set<Element> bag;
        for (Element e = 0; e < n; ++e) {
            // introduce elements
            for (auto ei: introduced_[e]) {
                bag.insert(ei);
            }
            // determine frontier
            for (auto eb: bag) {
                frontiers_[e].emplace_back(eb);
            }
            // forget elements
            for (auto ef: forgotten_[e]) {
                bag.erase(ef);
            }
        }
        assert(bag.empty());

        // preparare compared_
        compared_.resize(n);
        std::vector<bool> visit_compared(n, false);
        for (Element e = 0; e < n; ++e) {
            for (Element efr: frontiers_[e]){
                if (perm_.permute(efr) == e) {
                    assert(efr <= e);
                    // efr and e can be comapred.
                    compared_[e].emplace_back(efr);
                    visit_compared[efr] = true;
                }
            }
        }
        for (Element e = 0; e < n; ++e) {
            assert(visit_compared[e]);
        }

        // path with of the "propagation" graph
        max_frontier_size_ = 0;
        for (const auto& f: frontiers_) {
            max_frontier_size_ = std::max(static_cast<int>(f.size()), max_frontier_size_);
        }

        // map element to position
        mapping_element_to_pos_.resize(n);
        std::vector<bool> used(max_frontier_size_, false);
        for (Element e = 0; e < n; ++e) {
            // introduce elements
            for (auto ei: introduced_[e]) {
                bool success = false;
                for (FrontierPosition i = 0; i < max_frontier_size_; ++i) {
                    if (!used[i]) {
                        mapping_element_to_pos_[ei] = i;
                        used[i] = true;
                        success = true;
                        break;
                    }
                }
                assert(success);
            }
            // forget elements
            for (auto ef: forgotten_[e]) {
                FrontierPosition released = mapping_element_to_pos_[ef];
                used[released] = false;
            }
        }
        for (FrontierPosition i = 0; i < max_frontier_size_; ++i) {
            assert(!used[i]);
        }
    }
};

} // namespace permutation

#endif // PYZDD_PERMUTATION_H
