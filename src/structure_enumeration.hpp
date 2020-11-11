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
#include "spec/combination.hpp"
#include "spec/isomorphism.hpp"
#include "spec/superperiodic.hpp"
#include "spec/induced_subgraph.hpp"

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
        const std::vector<permutation::Permutation>& translations,
        const std::vector<int>& composition_constraints,
        const std::vector<std::vector<permutation::Element>>& site_constraints,
        bool remove_incomplete,
        bool remove_superperiodic,
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

        // spec for superperiodic elimination
        std::vector<permutation::superperiodic::SuperperiodicElimination> sp_specs;
        if (remove_superperiodic) {
            sp_specs.reserve(translations.size());
            auto identity = permutation::get_identity(num_variables);
            for (const auto& perm: translations) {
                if (perm == identity) {
                    continue;
                }
                permutation::PermutationFrontierManager pfm(perm);
                permutation::superperiodic::SuperperiodicElimination spec(pfm);
                sp_specs.emplace_back(spec);
            }
        }

        // spec for composition constraints
        std::vector<combination::Combination> composition_specs;
        if (!composition_constraints.empty()) {
            int k = composition_constraints[1]; // take k elements of 1-branchs
            combination::Combination spec(num_variables, k);
            composition_specs.emplace_back(spec);
        }

        // spec for site constraints
        std::vector<choice::Choice> site_constraint_specs;
        if (!site_constraints.empty()) {
            site_constraint_specs.reserve(num_sites);
            for (int site = 0; site < num_sites; ++site) {
                std::vector<int> group = {site};
                for (auto specie: site_constraints[site]) {
                    if (specie == 1) {
                        // cannot take 1-branch for the site
                        choice::Choice spec(num_variables, 0, group, false);
                        site_constraint_specs.emplace_back(spec);
                    } else {
                        // cannot take 0-branch for the site
                        choice::Choice spec(num_variables, 1, group, false);
                        site_constraint_specs.emplace_back(spec);
                    }
                }
            }
        }

        // spec for removing incompletes
        choice::TakeBoth incomplete_spec(num_variables);

        // ==== construct DD ====
        dd = universe::Universe(num_variables);

        if (remove_incomplete) {
            dd.zddSubset(incomplete_spec);
            dd.zddReduce();
        }

        if (!composition_constraints.empty()) {
            for (const auto& spec: composition_specs) {
                dd.zddSubset(spec);
                dd.zddReduce();
            }
        }

        if (!site_constraints.empty()) {
            for (const auto& spec: site_constraint_specs) {
                dd.zddSubset(spec);
                dd.zddReduce();
            }
        }

        if (remove_superperiodic) {
            for (const auto& spec: sp_specs) {
                dd.zddSubset(spec);
                dd.zddReduce();
            }
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
        const std::vector<permutation::Permutation>& translations,
        const std::vector<int>& composition_constraints,
        const std::vector<std::vector<permutation::Element>>& site_constraints,
        bool remove_incomplete,
        bool remove_superperiodic,
        tdzdd::DdStructure<2>& dd
    ) {
        // ==== one-hot encoding for num_types >= 3 ====
        int num_variables = num_sites * num_types;
        auto converter = SiteSpecieConverter(num_sites, num_types);

        std::vector<permutation::Permutation> automorphism_enc;
        automorphism_enc.reserve(automorphism.size());
        for (const auto& perm: automorphism) {
            automorphism_enc.emplace_back(converter.augment_permutation(perm));
        }

        std::vector<permutation::Permutation> translations_enc;
        translations_enc.reserve(translations_enc.size());
        for (const auto& perm: translations) {
            translations_enc.emplace_back(converter.augment_permutation(perm));
        }

        std::vector<int> site_constraints_enc;
        if (!site_constraints.empty()) {
            for (int site = 0; site < num_sites; ++site) {
                for (int specie: site_constraints[site]) {
                    auto prohibited = converter.encode(site, specie);
                    site_constraints_enc.emplace_back(prohibited);
                }
            }
        }

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

        // spec for superperiodic elimination
        std::vector<permutation::superperiodic::SuperperiodicElimination> sp_specs;
        if (remove_superperiodic) {
            sp_specs.reserve(translations.size());
            auto identity = permutation::get_identity(num_variables);
            for (const auto& perm: translations_enc) {
                if (perm == identity) {
                    continue;
                }
                permutation::PermutationFrontierManager pfm(perm);
                permutation::superperiodic::SuperperiodicElimination spec(pfm);
                sp_specs.emplace_back(spec);
            }
        }

        // spec for composition constraints
        std::vector<choice::Choice> composition_specs;
        if (!composition_constraints.empty()) {
            for (int specie = 0; specie < num_types; ++specie) {
                int k = composition_constraints[specie];
                std::vector<int> sites_with_same_specie;
                sites_with_same_specie.reserve(num_sites);
                for (int site = 0; site < num_sites; ++site) {
                    sites_with_same_specie.emplace_back(converter.encode(site, specie));
                }

                choice::Choice spec(num_variables, k, sites_with_same_specie, false);
                composition_specs.emplace_back(spec);
            }
        }

        // spec for site constraints
        std::vector<choice::Choice> site_constraint_specs;
        if (!site_constraints.empty()) {
            for (int prohibited: site_constraints_enc) {
                std::vector<int> group = {prohibited};
                // cannot choose `prohibited`
                choice::Choice spec(num_variables, 0, group, false);
                site_constraint_specs.emplace_back(spec);
            }
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

        if (!composition_constraints.empty()) {
            for (const auto& spec: composition_specs) {
                dd.zddSubset(spec);
                dd.zddReduce();
            }
        }

        if (!site_constraints.empty()) {
            for (const auto& spec: site_constraint_specs) {
                dd.zddSubset(spec);
                dd.zddReduce();
            }
        }

        if (remove_superperiodic) {
            for (const auto& spec: sp_specs) {
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
    ///            is a desired number of type-k.
    /// @param[in] site_constraints (Optional) site_constraints[i] is a list of
    ///            species prohibited to locate at site-i.
    /// @param[in] remove_superperiodic iff true, remove superperiodic structures
    /// @param[in] remove_incomplete iff true, remove incomplete structures,
    ///            whose kinds of species is less than num_types.
    /// @param[out] dd DD for output
    void construct_derivative_structures(
        tdzdd::DdStructure<2>& dd,
        int num_sites,
        int num_types,
        const std::vector<permutation::Permutation>& automorphism,
        const std::vector<permutation::Permutation>& translations,
        const std::vector<int>& composition_constraints,
        const std::vector<std::vector<permutation::Element>>& site_constraints,
        bool remove_incomplete,
        bool remove_superperiodic
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
        for (const auto& perm: translations) {
            if (perm.get_size() != static_cast<size_t>(num_sites)) {
                std::cerr << "The number of elements of permutation should be num_sites." << std::endl;
                exit(1);
            }
        }
        if (!composition_constraints.empty()) {
            if (composition_constraints.size() != static_cast<size_t>(num_types)) {
                std::cerr << "The size of composition_constraints should be num_types." << std::endl;
                exit(1);
            }
            if (std::accumulate(composition_constraints.begin(), composition_constraints.end(), 0) != num_sites) {
                std::cerr << "The sum of composition_consraints should be num_sites." << std::endl;
                exit(1);
            }
        }
        if (!site_constraints.empty()) {
            if (site_constraints.size() != static_cast<size_t>(num_sites)) {
                std::cerr << "The size of site_constraints should be num_sites." << std::endl;
                exit(1);
            }
            for (const auto& types: site_constraints) {
                for (auto k: types) {
                    if ((k < 0) || (k >= static_cast<permutation::Element>(num_types))) {
                        std::cerr << "Invalid specie: " << k << std::endl;
                        exit(1);
                    }
                }
                if (types.size() >= static_cast<size_t>(num_types)) {
                    std::cerr << "Impossible to satisfy site constraints!" << std::endl;
                    exit(1);
                }
            }
        }
        if (remove_superperiodic && translations.empty()) {
            std::cerr << "Translational group is required for removing superperiodic structures.";
            exit(1);
        }

        if (num_types == 2) {
            enumerate_binary_derivative_structures(
                num_sites,
                num_types,
                automorphism,
                translations,
                composition_constraints,
                site_constraints,
                remove_incomplete,
                remove_superperiodic,
                dd
            );
        } else {
            enumerate_multi_derivative_structures(
                num_sites,
                num_types,
                automorphism,
                translations,
                composition_constraints,
                site_constraints,
                remove_incomplete,
                remove_superperiodic,
                dd
            );
        }
    }

    /// @brief enumerate DD for derivative structures with fixed Warren-Cowley SRO
    /// @param[out] dd DD for output
    /// @param[in] num_sites the number of sites in supercell
    /// @param[in] num_types the kinds of species
    /// @param[in] automorphism symmetry group of supercell
    /// @param[in] translations (Optional) permutation group derived from
    ////           translations of cell. Required if `remove_superperiodic` is true.
    /// @param[in] composition_constraints composition_constraints[i]
    ///            is a pair of sites and a desired number of label=1
    /// @param[in] vgfm frontier manager of cluster graph
    /// @param[in] target target weight without loop offset of vertex-induced subgraph
    void construct_binary_derivative_structures_with_sro(
        tdzdd::DdStructure<2>& dd,
        int num_sites,
        int num_types,
        const std::vector<permutation::Permutation>& automorphism,
        const std::vector<permutation::Permutation>& translations,
        const std::vector<std::pair<std::vector<int>, int>>& composition_constraints,
        const graph::VertexGraphFrontierManager& vgfm,
        graph::Weight target
    ) {
        assert(num_types == 2);
        size_t num_variables = num_sites;

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

        // spec for superperiodic elimination
        std::vector<permutation::superperiodic::SuperperiodicElimination> sp_specs;
        sp_specs.reserve(translations.size());
        auto identity = permutation::get_identity(num_variables);
        for (const auto& perm: translations) {
            if (perm == identity) {
                continue;
            }
            permutation::PermutationFrontierManager pfm(perm);
            permutation::superperiodic::SuperperiodicElimination spec(pfm);
            sp_specs.emplace_back(spec);
        }

        // spec for composition constraints
        assert(!composition_constraints.empty());
        std::vector<choice::Choice> composition_specs;
        composition_specs.reserve(composition_constraints.size());
        for (const auto& variables_count: composition_constraints) {
            std::vector<int> group;
            group.reserve(variables_count.first.size());
            for (graph::Vertex v: variables_count.first) {
                group.emplace_back(vgfm.map_to_internal_vertex_id(v));
            }

            int k = variables_count.second; // take k elements of 1-branchs
            choice::Choice spec(num_variables, k, group, false);
            composition_specs.emplace_back(spec);
        }

        // spec for SRO
        graph::induced_subgraph::VertexInducedSubgraphSpec subgraph_spec(vgfm, target);

        // spec for removing incompletes
        choice::TakeBoth incomplete_spec(num_variables);

        // ==== construct DD ====
        dd = universe::Universe(num_variables);

        // remove incomplete
        dd.zddSubset(incomplete_spec);
        dd.zddReduce();

        // fix composition
        for (const auto& spec: composition_specs) {
            dd.zddSubset(spec);
            dd.zddReduce();
        }

        // remove superperiodic
        for (const auto& spec: sp_specs) {
            dd.zddSubset(spec);
            dd.zddReduce();
        }

        // remove symmetry duplicates
        for (const auto& spec: aut_specs) {
            dd.zddSubset(spec);
            dd.zddReduce();
        }

        // fix SRO
        dd.zddSubset(subgraph_spec);
        dd.zddReduce();
    }

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

    std::vector<int> convert_to_labeling_with_graph(
        const tdzdd::DdStructure<2>::const_iterator& itr,
        const graph::VertexGraphFrontierManager& vgfm,
        int num_types
    )
    {
        assert(num_types == 2);
        // vertex order in DD is diffenrent from the original variable order in the graph
        std::vector<int> labeling = vgfm.retrieve_vertices(*itr);
        return labeling;
    }


} // namespace derivative_structure
} // namespace pyzdd

#endif // PYZDD_STRUCTURE_ENUMERATION_H
