#ifndef PYZDD_SHORT_RANGE_ORDRE_H
#define PYZDD_SHORT_RANGE_ORDRE_H

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
    const std::vector<graph::VertexGraphFrontierManager>& vgfm_vec,
    const std::vector<graph::Weight>& targets
) {
    tdzdd::MessageHandler::showMessages(true);

    assert(num_types == 2);
    assert(!vgfm_vec.empty());
    size_t num_variables = num_sites;

    // ==== translate permutations with vertex_order ====
    const auto& vgfm0 = vgfm_vec[0];

    std::vector<permutation::Permutation> reordered_automorphism;
    reordered_automorphism.reserve(automorphism.size());
    for (const auto& perm: automorphism) {
        std::vector<permutation::Element> sigma(num_variables);
        for (graph::Vertex v = 0; v < static_cast<graph::Vertex>(num_variables); ++v) {
            sigma[vgfm0.map_to_internal_vertex_id(v)] = vgfm0.map_to_internal_vertex_id(perm.permute(v));
        }
        reordered_automorphism.emplace_back(permutation::Permutation(sigma));
    }
    // sort permutations by max frontier sizes
    std::sort(reordered_automorphism.begin(), reordered_automorphism.end(),
                [](const permutation::Permutation& lhs, const permutation::Permutation& rhs) {
                    auto size_lhs = permutation::PermutationFrontierManager(lhs).get_max_frontier_size();
                    auto size_rhs = permutation::PermutationFrontierManager(rhs).get_max_frontier_size();
                    return size_lhs < size_rhs;
                });

    std::vector<permutation::Permutation> reordered_translation_group;
    reordered_translation_group.reserve(translations.size());
    for (const auto& perm: translations) {
        std::vector<permutation::Element> sigma(num_variables);
        for (graph::Vertex v = 0; v < static_cast<graph::Vertex>(num_variables); ++v) {
            sigma[vgfm0.map_to_internal_vertex_id(v)] = vgfm0.map_to_internal_vertex_id(perm.permute(v));
        }
        reordered_translation_group.emplace_back(permutation::Permutation(sigma));
    }

    // ==== prepare specs ====
    // spec for isomorphism elimination
    std::vector<permutation::Permutation> sorted_automorphism(automorphism);

    std::vector<permutation::isomorphism::IsomorphismElimination> aut_specs;
    aut_specs.reserve(reordered_automorphism.size());
    for (const auto& perm: reordered_automorphism) {
        permutation::PermutationFrontierManager pfm(perm);
        permutation::isomorphism::IsomorphismElimination spec(pfm);
        aut_specs.emplace_back(spec);
    }

    // spec for superperiodic elimination
    std::vector<permutation::superperiodic::SuperperiodicElimination> sp_specs;
    sp_specs.reserve(reordered_translation_group.size());
    auto identity = permutation::get_identity(num_variables);
    for (const auto& perm: reordered_translation_group) {
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
            group.emplace_back(vgfm0.map_to_internal_vertex_id(v));
        }

        int k = variables_count.second; // take k elements of 1-branchs
        choice::Choice spec(num_variables, k, group, false);
        composition_specs.emplace_back(spec);
    }

    // spec for SRO
    std::vector<graph::induced_subgraph::VertexInducedSubgraphSpec> subgraph_specs;
    size_t num_graphs = vgfm_vec.size();
    assert(targets.size() == num_graphs);
    if (num_graphs > 0) {
        subgraph_specs.reserve(num_graphs);
        // TODO: better vertex_order choice
        const auto& vertex_order = vgfm_vec[0].get_vertex_order();

        for (size_t i = 0; i < num_graphs; ++i) {
            // need to use the same vertex-order!
            assert(vgfm_vec[i].get_vertex_order() == vertex_order);

            graph::induced_subgraph::VertexInducedSubgraphSpec spec(vgfm_vec[i], targets[i]);
            subgraph_specs.emplace_back(spec);
        }
    }

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

    // fix SRO before uniquing by symmetry
    for (const auto& spec: subgraph_specs) {
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

    tdzdd::MessageHandler::showMessages(false);
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

#endif // PYZDD_SHORT_RANGE_ORDRE_H
