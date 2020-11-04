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

    class SiteSpecieConverter {
        const size_t num_sites_;
        const size_t num_types_;
        const size_t num_variables_;
    public:
        SiteSpecieConverter() = delete;
        SiteSpecieConverter(int num_sites, int num_types) :
            num_sites_(static_cast<size_t>(num_sites)),
            num_types_(static_cast<size_t>(num_types)),
            num_variables_(num_sites * num_types)
        {
            assert(num_types_ >= 2);
        }

        /// @brief encode (site, specie)
        permutation::Element encode(permutation::Element site, int specie) const {
            return site * num_types_ + specie;
        }

        /// @brief return (site, specie)
        std::pair<permutation::Element, int> decode(permutation::Element idx) const {
            return std::make_pair(idx / num_types_, idx % num_types_);
        }

        permutation::Permutation augment_permutation(const permutation::Permutation& perm) const {
            std::vector<permutation::Element> sigma(num_variables_);
            for (size_t site = 0; site < num_sites_; ++site) {
                for (size_t specie = 0; specie < num_types_; ++specie) {
                    sigma[encode(perm.permute(site), specie)] = encode(site, specie);
                }
            }
            auto ret = permutation::Permutation(sigma);
            return ret;
        }
    };

    /// @brief enumerate binary (num_types == 2) derivative structures
    void enumerate_binary_derivative_structures(
        int num_sites,
        int num_types,
        const std::vector<permutation::Permutation>& automorphism,
        bool remove_incomplete,
        tdzdd::DdStructure<2>& dd
    ) {
        int num_variables = num_sites;

        // ==== prepare specs ====
        // spec for isomorphism elimination
        // sort permutations by max frontier sizes
        std::vector<permutation::Permutation> sorted_automorphism(automorphism);
        std::sort(sorted_automorphism.begin(), sorted_automorphism.end(),
                  [](const permutation::Permutation& lhs, const permutation::Permutation& rhs) {
                      auto size_lhs = permutation::PermutationFrontierManager(lhs).get_max_frontier_size();
                      auto size_rhs = permutation::PermutationFrontierManager(rhs).get_max_frontier_size();
                      return size_lhs < size_rhs;
                  });
        std::vector<permutation::isomorphism::IsomorphismElimination> aut_specs;
        aut_specs.reserve(sorted_automorphism.size());
        for (const auto& perm: sorted_automorphism) {
            permutation::PermutationFrontierManager pfm(perm);
            permutation::isomorphism::IsomorphismElimination spec(pfm);
            aut_specs.emplace_back(spec);
        }

        // spec for removing incompletes
        choice::TakeBoth incomplete_spec(num_variables);

        // ==== construct DD ====
        dd = universe::Universe(num_variables);

        if (remove_incomplete) {
            dd.zddSubset(incomplete_spec);
            dd.zddReduce();
        }

        for (const auto& spec: aut_specs) {
            dd.zddSubset(spec);
            dd.zddReduce();
        }
    }

    /// @brief enumerate binary (num_types == 2) derivative structures
    void enumerate_multi_derivative_structures(
        int num_sites,
        int num_types,
        const std::vector<permutation::Permutation>& automorphism,
        bool remove_incomplete,
        tdzdd::DdStructure<2>& dd
    ) {
        // one-hot encoding for num_types >= 3
        int num_variables = num_sites * num_types;
        auto converter = SiteSpecieConverter(num_sites, num_types);

        std::vector<permutation::Permutation> automorphism_enc;
        automorphism_enc.reserve(automorphism.size());
        for (const auto& perm: automorphism) {
            automorphism_enc.emplace_back(converter.augment_permutation(perm));
        }

        // std::vector<permutation::Permutation> translations_enc;
        // std::vector<std::pair<std::vector<int>, int>> composition_constraints_enc;

        // ==== prepare specs ====
        // spec for isomorphism elimination
        // sort permutations by max frontier sizes
        std::sort(automorphism_enc.begin(), automorphism_enc.end(),
                  [](const permutation::Permutation& lhs, const permutation::Permutation& rhs) {
                      auto size_lhs = permutation::PermutationFrontierManager(lhs).get_max_frontier_size();
                      auto size_rhs = permutation::PermutationFrontierManager(rhs).get_max_frontier_size();
                      return size_lhs < size_rhs;
                  });
        std::vector<permutation::isomorphism::IsomorphismElimination> aut_specs;
        aut_specs.reserve(automorphism_enc.size());
        for (const auto& perm: automorphism_enc) {
            permutation::PermutationFrontierManager pfm(perm);
            permutation::isomorphism::IsomorphismElimination spec(pfm);
            aut_specs.emplace_back(spec);
        }

        // spec for removing incompletes
        std::vector<choice::Choice> incomplete_specs;
        incomplete_specs.reserve(num_types);
        for (int specie = 0; specie < num_types; ++specie) {
            // elements with the same specie
            std::vector<int> group;
            group.reserve(num_sites);
            for (int site = 0; site < num_sites; ++site) {
                group.emplace_back(converter.encode(site, specie));
            }

            // take elements in `group` at least 1
            choice::Choice spec(num_variables, 1, group, true);
            incomplete_specs.emplace_back(spec);
        }

        // spec for one-of-k
        std::vector<choice::Choice> onehot_specs;
        onehot_specs.reserve(num_sites);
        for (permutation::Element site = 0; site < static_cast<permutation::Element>(num_sites); ++site) {
            std::vector<int> same_sites(num_types);
            for (int specie = 0; specie < num_types; ++specie) {
                same_sites[specie] = converter.encode(site, specie);
            }
            choice::Choice spec(num_variables, 1, same_sites, false);
            onehot_specs.emplace_back(spec);
        }

        // ==== construct DD ====
        dd = universe::Universe(num_variables);
        for (const auto& spec: onehot_specs) {
            dd.zddSubset(spec);
            dd.zddReduce();
        }

        if (remove_incomplete) {
            for (const auto& spec: incomplete_specs) {
                dd.zddSubset(spec);
                dd.zddReduce();
            }
        }

        for (const auto& spec: aut_specs) {
            dd.zddSubset(spec);
            dd.zddReduce();
        }
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
        bool remove_incomplete,
        tdzdd::DdStructure<2>& dd
    ) {
        // sanity check
        assert(num_sites >= 1);
        assert(num_types >= 2);
        for (const auto& perm: automorphism) {
            if (perm.get_size() != static_cast<size_t>(num_sites)) {
                std::cerr << "The number of elements of permutation should be num_sites." << std::endl;
                exit(1);
            }
        }

        if (num_types == 2) {
            enumerate_binary_derivative_structures(num_sites, num_types, automorphism, remove_incomplete, dd);
        } else {
            enumerate_multi_derivative_structures(num_sites, num_types, automorphism, remove_incomplete, dd);
        }
    }

    /*
    void enumerate_derivative_structures(
        const std::vector<permutation::Permutation>& translations,
        const std::vector<std::pair<std::vector<int>, int>>& composition_constraints,
        const std::vector<std::vector<permutation::Element>>& site_constraints,
        bool remove_superperiodic,
    ) {
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


        // std::vector<std::vector<permutation::Element>> prohibited_species_enc;
        if (num_types == 2) {
            if (!site_constraints.empty()) {
                prohibited_species_enc.reserve(site_constraints.size());
                for (permutation::Element site = 0; site < num_sites; ++site) {
                    std::set<int> allow;
                    for (auto specie: site_constraints[site]) {
                        allow.insert(specie);
                    }

                    std::vector<int> prohibited;
                    for (int specie = 0; specie < num_types; ++specie) {
                        if (allow.find(specie) == allow.end()) {
                            prohibited.emplace_back(specie);
                        }
                    }
                    prohibited_species_enc.emplace_back(prohibited);
                }
            }
        } else {

            translations_enc.reserve(translations.size());
            for (const auto& perm: translations) {
                translations_enc.emplace_back(converter.augment_permutation(perm));
            }

            if (!composition_constraints.empty()) {
                composition_constraints_enc.reserve(composition_constraints.size());
                for (int specie = 0; specie < num_types; ++specie) {
                    const auto& group = composition_constraints[specie].first;
                    int target_number = composition_constraints[specie].second;
                    std::vector<int> augmented_group;
                    augmented_group.reserve(group.size());
                    for (auto site: group) {
                        augmented_group.emplace_back(converter.encode(site, specie));
                    }
                    composition_constraints_enc.emplace_back(std::make_pair(augmented_group, target_number));
                }
            }

            if (!site_constraints.empty()) {
                prohibited_species_enc.reserve(site_constraints.size());
                for (permutation::Element site = 0; site < num_sites; ++site) {
                    std::set<int> allow;
                    for (auto specie: site_constraints[site]) {
                        allow.insert(specie);
                    }

                    std::vector<int> prohibited;
                    for (int specie = 0; specie < num_types; ++specie) {
                        if (allow.find(specie) == allow.end()) {
                            prohibited.emplace_back(specie);
                        }
                    }

                    std::vector<permutation::Element> augmented;
                    augmented.reserve(site_constraints[site].size());
                    for (auto specie: site_constraints[site]) {
                        augmented.emplace_back(converter.encode(site, specie));
                    }
                    site_constraints_enc.emplace_back(augmented);
                }
            }
        }


        // spec for superperiodic elimination
        std::vector<permutation::superperiodic::SuperperiodicElimination> sp_specs;
        sp_specs.reserve(translations_enc.size());
        for (const auto& perm: translations_enc) {
            permutation::PermutationFrontierManager pfm(perm);
            permutation::superperiodic::SuperperiodicElimination spec(pfm);
            sp_specs.emplace_back(spec);
        }

        // spec for composition constraints
        std::vector<choice::Choice> composition_specs;
        if (!composition_constraints_enc.empty()) {
            composition_specs.reserve(num_types);
            for (const auto& pair: composition_constraints_enc) {
                const auto& group = pair.first;
                int k = pair.second;
                // take exactly k elements from group (do not care element not in group)
                choice::Choice spec(num_variables, k, group, false);
                composition_specs.emplace_back(spec);
            }
        }

        // spec for site constraints
        std::vector<pyzdd::choice::Choice> sites_specs;
        if (!site_constraints_enc.empty()) {
            sites_specs.reserve(num_variables);
            for (permutation::Element site = 0; site < num_sites; ++site) {
                std::set<int> allow;
                for (auto e: site_constraints_enc[site]) {

                }
            }

            for (const auto& group: site_constraints_enc) {
                // take exactly k elements from group (do not care element not in group)
                pyzdd::choice::Choice spec(num_variables, k, group, false);
                choice_specs.emplace_back(spec);
            }
        }

        // spec for one-hot encoding
    }
    */

    std::vector<permutation::Element> convert_to_labeling(
        tdzdd::DdStructure<2>::const_iterator const &itr,
        int num_sites,
        int num_types)
    {
        int num_variables = (num_types == 2) ? num_sites : (num_sites * num_types);
        std::vector<permutation::Element> labeling(num_sites, 0);

        if (num_types == 2) {
            for (auto level: *itr) {
                labeling[num_variables - level] = 1;
            }
        } else {
            auto converter = SiteSpecieConverter(num_sites, num_types);
            for (auto level: *itr) {
                auto site_specie = converter.decode(num_variables - level);
                labeling[site_specie.first] = site_specie.second;
            }
        }

        return labeling;
    }

} // namespace derivative_structure
} // namespace pyzdd

#endif // PYZDD_STRUCTURE_ENUMERATION_H
