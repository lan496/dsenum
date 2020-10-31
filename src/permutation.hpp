#ifndef PYZDD_PERMUTATION_H
#define PYZDD_PERMUTATION_H

#include <iostream>
#include <algorithm>
#include <vector>
#include <set>
#include <queue>
#include <numeric>
#include <utility>
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

/// @brief Permutaion class
class Permutation {
    /// the number of elements
    Element n_;
    /// represents a permutation in "one-line" notation. That is, `sigma` moves `i` to `sigma[i]`.
    std::vector<Element> sigma_;
public:
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
            std::cerr << "Given list is not a permutation:";
            for (auto e: sigma) {
                std::cerr << " " << static_cast<int>(e);
            }
            std::cerr << std::endl;
            exit(1);
        }
    }

    Permutation() = default;
    Permutation(const Permutation&) = default;
    Permutation& operator=(const Permutation&) = default;
    Permutation(Permutation&& other) = default;
    Permutation& operator=(Permutation&&) = default;

    /// @brief get the number of elements
    size_t get_size() const {
        return n_;
    }

    /// @brief permute i to sigma[i]
    Element permute(const Element i) const {
        return sigma_[i];
    }

    /// @brief act the permutation on a given sequence
    /// @tparam T is template param
    /// @param[in] colors colors[i] is permutated to colors[perm-1[i]]
    template<typename T>
    std::vector<T> act(const std::vector<T> colors) const {
        std::vector<T> permutated(colors.size());
        for (Element i = 0, n = get_size(); i < n; ++i) {
            permutated[i] = colors[permute(i)];
        }
        return permutated;
    }

    /// @brief self.product(rhs).permute(i) == self.permute(rhs.permute(i))
    Permutation product(const Permutation& rhs) const {
        auto n = get_size();
        if (rhs.get_size() != n) {
            std::cerr << "Cannot product permutations with the different bases." << std::endl;
            exit(1);
        }
        std::vector<Element> sigma(n);
        for (Element i = 0; i < n; ++i) {
            sigma[i] = permute(rhs.permute(i));
        }
        return Permutation(sigma);
    }

    /// @brief return inverse of the permutation
    Permutation inverse() const {
        auto n = get_size();
        std::vector<Element> sigma(n);
        for (Element i = 0; i < n; ++i) {
            sigma[permute(i)] = i;
        }
        return Permutation(sigma);
    }

    void dump(std::ostream& os) const {
        os << "(";
        for (Element i = 0, n = get_size(); i < n; ++i) {
            os << " " << permute(i);
        }
        os << " )" << std::endl;
    }

private:
    /// @brief check if `sigma` is a bijection.
    bool is_permutation(const std::vector<Element>& sigma) const {
        auto n = get_size();
        std::vector<bool> visited(n, false);
        for (size_t i = 0; i < n; ++i) {
            if ((sigma[i] >= n) || visited[sigma[i]]) {
                return false;
            }
            visited[sigma[i]] = true;
        }
        return true;
    }
};

inline bool operator==(const Permutation& lhs, const Permutation& rhs) {
    if (lhs.get_size() != rhs.get_size()) {
        return false;
    }
    for (Element i = 0, n = lhs.get_size(); i < n; ++i) {
        if (lhs.permute(i) != rhs.permute(i)) {
            return false;
        }
    }
    return true;
}

inline bool operator!=(const Permutation& lhs, const Permutation& rhs) {
    return !(lhs == rhs);
}

/// @brief return an identity permutation on `n` elements
Permutation get_identity(Element n) {
    std::vector<Element> sigma(n);
    for (Element i = 0; i < n; ++i) {
        sigma[i] = i;
    }
    return Permutation(sigma);
}

/// @brief generate a permutation group from given generators
std::vector<Permutation> generate_group(const std::vector<Permutation>& generators) {
    if (generators.empty()) {
        return std::vector<Permutation>();
    }

    auto n = generators[0].get_size();
    for (auto g: generators) {
        assert(g.get_size() == n);
    }

    std::vector<Permutation> group;
    std::queue<Permutation> que;
    auto id = get_identity(n);
    group.emplace_back(id);
    que.push(id);
    while (!que.empty()) {
        auto q = que.front(); que.pop();
        for (auto g: generators) {
            // from left
            auto gq = g.product(q);
            if (std::find(group.begin(), group.end(), gq) == group.end()) {
                group.emplace_back(gq);
                que.push(gq);
            }
            // from right
            auto qg = q.product(g);
            if (std::find(group.begin(), group.end(), qg) == group.end()) {
                group.emplace_back(qg);
                que.push(qg);
            }
        }
    }
    return group;
}

// ============================================================================
// Prepare frontiers for Isomorphism Elimination
// ============================================================================

/// @brief construct and manage frontiers of a permutation
class PermutationFrontierManager {
    /// permutation
    Permutation perm_;
    /// inverse of the permutation
    Permutation inv_perm_;

    /// Just before coloring the i-th element, the states of elements in
    /// frontiers_[i] are required. Guaranteed to be sorted in ascending order.
    std::vector<std::vector<Element>> frontiers_;
    /// Just after coloring the i-th element, original element compared_[i][].first
    /// and permuted element compared_[i][].second can be compared.
    std::vector<std::vector<std::pair<Element,Element>>> compared_;
    /// After processing the i-th element, color of elements in forgotten_[i] are no more required.
    std::vector<std::vector<Element>> forgotten_;
    /// mapping_element_to_pos_[e] is a position of element e in frontiers.
    /// Be careful several elements may share the same position!
    std::vector<FrontierPosition> mapping_element_to_pos_;
    /// the maximum size of frontiers
    int max_frontier_size_;

    /// After processing the i-th element, for j in comp_finished_[i], the
    /// comparison with element j and sigma(j) is already finished.
    /// comp_finished_[i] is sorted in the ascending order.
    std::vector<std::vector<Element>> comp_finished_;

public:
    PermutationFrontierManager() = delete;
    PermutationFrontierManager(const Permutation& perm) :
        perm_(perm),
        inv_perm_(perm.inverse())
    {
        // construct frontiers, forgotten, and max_frontier_size
        construct_coloring_frontier();
        construct_comparison_frontier();
    }

    /// @brief get the number of elements of the permutation
    Element get_size() const {
        return perm_.get_size();
    }

    int get_max_frontier_size() const {
        return max_frontier_size_;
    }

    const std::vector<Element>& get_frontier(Element e) const {
        return frontiers_[e];
    }

    const std::vector<std::pair<Element, Element>>& get_compared(Element e) const {
        return compared_[e];
    }

    const std::vector<Element>& get_forgotten(Element e) const {
        return forgotten_[e];
    }

    const FrontierPosition map_element_to_position(Element e) const {
        return mapping_element_to_pos_[e];
    }

    const std::vector<Element>& get_comp_finished(Element e) const {
        return comp_finished_[e];
    }

    /// @brief dump all frontier-search related information
    void dump(std::ostream& os) const {
        // show permutation
        os << "permutation" << std::endl;
        perm_.dump(os);

        // show frontier related
        auto n = get_size();

        os << "frontiers" << std::endl;
        for (Element e = 0; e < n; ++e) {
            os << "     " << e << ":";
            for (auto ef: get_frontier(e)) {
                os << " " << ef;
            }
            os << std::endl;
        }

        os << "compared" << std::endl;
        for (Element e = 0; e < n; ++e) {
            os << "     " << e << ":";
            for (auto pair: get_compared(e)) {
                os << " (" << pair.first << ", " << pair.second << ")";
            }
            os << std::endl;
        }

        os << "forgotton" << std::endl;
        for (Element e = 0; e < n; ++e) {
            os << "     " << e << ":";
            for (auto ef: get_forgotten(e)) {
                os << " " << ef;
            }
            os << std::endl;
        }

        os << "mapping" << std::endl;
        for (Element e = 0; e < n; ++e) {
            os << " " << map_element_to_position(e);
        }
        os << std::endl;

        os << "finished comparison" << std::endl;
        for (Element e = 0; e < n; ++e) {
            os << "     " << e << ":";
            for (auto efc: get_comp_finished(e)) {
                os << " " << efc;
            }
            os << std::endl;
        }

        os << std::endl;
    }
private:
    /// @brief construct frontiers and related
    void construct_coloring_frontier() {
        auto n = get_size();

        // prepare forgotten_
        // We can forget a Element e such that both (e, sigma(e)) and
        // (sigma-1(e), e) are already compared.
        forgotten_.resize(n);
        for (Element e = 0; e < n; ++e) {
            auto lifetime = std::max(e, std::max(perm_.permute(e), inv_perm_.permute(e)));
            forgotten_[lifetime].emplace_back(e);
        }

        // prepare frontiers
        frontiers_.resize(n);
        std::set<Element> bag;
        for (Element e = 0; e < n; ++e) {
            // introduce elements
            bag.insert(e);

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

        // prepare compared_
        // In comparison phase, (e, sigma(e)) and (simga-1(e), e) could be compared.
        compared_.resize(n);
        std::vector<bool> visit_compared(n, false);
        for (Element e = 0; e < n; ++e) {
            auto e_perm = perm_.permute(e);
            if (e_perm == e) {
                // e is isolated
                compared_[e].emplace_back(std::make_pair(e, e));
                assert(!visit_compared[e]);
                visit_compared[e] = true;
            } else {
                // compare e and sigma(e)
                if (e_perm < e) {
                    compared_[e].emplace_back(std::make_pair(e, e_perm));
                    assert(!visit_compared[e]);
                    visit_compared[e] = true;
                }
                // compare sigma-1(e) and e
                auto e_invperm = inv_perm_.permute(e);
                if (e_invperm < e) {
                    compared_[e].emplace_back(std::make_pair(e_invperm, e));
                    // e_invperm should be contained in the frontier
                    assert(std::find(frontiers_[e].begin(), frontiers_[e].end(), e_invperm) != frontiers_[e].end());
                    assert(!visit_compared[e_invperm]);
                    visit_compared[e_invperm] = true;
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
            bool success = false;
            for (FrontierPosition i = 0; i < static_cast<FrontierPosition>(max_frontier_size_); ++i) {
                if (!used[i]) {
                    mapping_element_to_pos_[e] = i;
                    used[i] = true;
                    success = true;
                    break;
                }
            }
            assert(success);

            // forget elements
            for (auto ef: forgotten_[e]) {
                FrontierPosition released = mapping_element_to_pos_[ef];
                used[released] = false;
            }
        }
        for (FrontierPosition i = 0; i < static_cast<FrontierPosition>(max_frontier_size_); ++i) {
            assert(!used[i]);
        }
    }

    void construct_comparison_frontier() {
        auto n = get_size();
        comp_finished_.resize(n);
        std::set<Element> finished;

        for (Element e = 0; e < n; ++e) {
            const auto& compared = get_compared(e);
            for (auto pair: compared) {
                // comparison (ef, sigma(ef)) is executed after determining the value of element e
                Element ef = pair.first;
                finished.insert(ef);
            }
            comp_finished_[e].reserve(finished.size());
            // insert with sorting
            for (auto itr = finished.begin(), end =finished.end(); itr != end; ++itr) {
                comp_finished_[e].emplace_back(*itr);
            }
        }
    }
};

} // namespace permutation

#endif // PYZDD_PERMUTATION_H
