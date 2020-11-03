#ifndef PYZDD_STRUCTURE_ENUMERATION_H
#define PYZDD_STRUCTURE_ENUMERATION_H

#include <vector>
#include <utility>

#include <tdzdd/DdStructure.hpp>
#include "type.hpp"
#include "iterator.hpp"
#include "permutation.hpp"
#include "spec/universe.hpp"
#include "spec/choice.hpp"
#include "spec/isomorphism.hpp"
#include "spec/superperiodic.hpp"

namespace pyzdd {

namespace derivative_structure {
    /// @brief encode (site, specie)
    int encode_site_with_specie(int site, int specie, int num_types) {
        return site * num_types + specie;
    }

    /// @brief return (site, specie)
    std::pair<int, int> decode_to_site_and_specie(int idx, int num_types) {
        return std::make_pair(idx / num_types, idx % num_types);
    }

    std::vector<permutation::Permutation> augment_permutations(const std::vector<permutation::Permutation>& perms, int num_sites, int num_types) {
        std::vector<permutation::Permutation> augmented;
        if (num_types == 2) {
            // binary system
            augmented = perms;
        } else {
            augmented.reserve(perms.size());
            for (const auto& perm: perms) {
                std::vector<permutation::Element> sigma(num_sites * num_types);
                for (size_t site = 0; site < num_sites; ++site) {
                    for (size_t specie = 0; specie < num_types; ++specie) {
                        sigma[encode_site_with_specie(perm.permute(site), specie, num_types)] = encode_site_with_specie(site, specie, num_types);
                    }
                }
                augmented.emplace_back(permutation::Permutation(sigma));
            }
        }
        return augmented;
    }

    /// @brief enumerate DD for nonequivalent derivative structures
    /// @details see K. Shinohara, et al., J. Chem. Phys. 153, 104109 (2020).
    /// @param[in] num_sites the number of sites in supercell
    /// @param[in] num_types the kinds of species
    /// @param[in] automorphism symmetry group of supercell
    /// @param[in] translations (Optional) permutation group derived from
    ////           translations of cell. Required if `remove_superperiodic` is true.
    /// @param[in] composition_constraints (Optional) composition_constraints[k]
    ///            is a pair of sites for type-k and the desired number of type-k.
    /// @param[in] site_constraints (Optional) site_constraints[i] is a list of
    ///            species allowed to locate at site-i.
    /// @param[in] remove_superperiodic iff true, remove superperiodic structures
    /// @param[in] remove_incomplete iff true, remove incomplete structures,
    ///            whose kinds of species is less than num_types.
    /// @param[out] dd DD for output
    void enumerate_derivative_structures(
        int num_sites,
        int num_types,
        const std::vector<permutation::Permutation>& automorphism,
        const std::vector<permutation::Permutation>& translations,
        const std::vector<std::pair<std::vector<int>, int>>& composition_constraints,
        const std::vector<std::vector<int>>& site_constraints,
        bool remove_superperiodic,
        bool remove_incomplete,
        tdzdd::DdStructure<2>& dd
    ) {
        // sanity check
        assert(num_sites >= 1);
        assert(num_types >= 2);
        for (const auto& perm: automorphism) {
            if (perm.get_size() != num_sites) {
                std::cerr << "The number of elements of permutation should be num_sites." << std::endl;
                exit(1);
            }
        }
        for (const auto& perm: translations) {
            if (perm.get_size() != num_sites) {
                std::cerr << "The number of elements of permutation should be num_sites." << std::endl;
                exit(1);
            }
        }
        if (!composition_constraints.empty()) {
            if (composition_constraints.size() != static_cast<size_t>(num_types)) {
                std::cerr << "The size of composition_constraints should be num_types." << std::endl;
                exit(1);
            }
            for (const auto& pair: composition_constraints) {
                if (pair.second < 0) {
                    std::cerr << "The concentration should not be negative." << std::endl;
                    exit(1);
                }
            }
        }
        if (!site_constraints.empty()) {
            if (site_constraints.size() != static_cast<size_t>(num_sites)) {
                std::cerr << "The size of site_constraints should be num_sites." << std::endl;
                exit(1);
            }
            for (const auto& types: site_constraints) {
                for (auto k: types) {
                    if ((k < 0) || (k >= num_types)) {
                        std::cerr << "Invalid specie: " << k << std::endl;
                        exit(1);
                    }
                }
            }
        }
        if (remove_superperiodic && translations.empty()) {
            std::cerr << "Translational group is required for removing superperiodic structures.";
            exit(1);
        }

        // one-hot encoding for num_types >= 3
        bool is_binary_system = (num_types == 2);
        int num_variables = is_binary_system ? num_sites : (num_sites * num_types);

        auto automorphism_enc = augment_permutations(automorphism, num_sites, num_types);
        auto translations_enc = augment_permutations(translations, num_sites, num_types);

        std::vector<std::pair<std::vector<int>, int>> composition_constraints_enc;
        if (!composition_constraints.empty()) {
            if (is_binary_system) {
                composition_constraints_enc = composition_constraints;
            } else {
                composition_constraints_enc.reserve(composition_constraints.size());
                for (int specie = 0; specie < num_types; ++num_types) {
                    // TODO
                }
            }
        }

        // prepare specs
        std::vector<pyzdd::choice::Choice> choice_specs;
        if (!composition_constraints.empty()) {
            // TODO
            choice_specs.reserve(num_types);
            for (const auto& pair: composition_constraints) {
                const auto& group = pair.first;
                int k = pair.second;
                // take exactly k elements from group (do not care element not in group)
                pyzdd::choice::Choice spec(num_variables, k, group, false);
                choice_specs.emplace_back(spec);
            }
        }

        tdzdd::DdStructure<2> dd = universe::Universe(num_variables);
    }

} // namespace derivative_structure
} // namespace pyzdd

#endif // PYZDD_STRUCTURE_ENUMERATION_H
