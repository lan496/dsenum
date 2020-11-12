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

class BinaryVertexConverter {
    size_t num_variables_;
    std::vector<graph::Vertex> vertex_order_;
    std::vector<graph::InternalVertexId> mapping_vertex_;
public:
    BinaryVertexConverter() = delete;
    BinaryVertexConverter(int num_sites, const std::vector<graph::Vertex>& vertex_order) :
        num_variables_(static_cast<size_t>(num_sites)),
        vertex_order_(vertex_order)
    {
        // mapping_vertex_
        mapping_vertex_.resize(num_variables_);
        for (graph::InternalVertexId vid = 0; vid < static_cast<graph::InternalVertexId>(num_variables_); ++vid) {
            mapping_vertex_[vertex_order_[vid]] = vid;
        }
    }

    permutation::Permutation reorder_premutation(const permutation::Permutation& perm) const {
        std::vector<permutation::Element> sigma(num_variables_);
        for (graph::Vertex v = 0; v < static_cast<graph::Vertex>(num_variables_); ++v) {
            sigma[map_to_internal_vertex_id(v)] = map_to_internal_vertex_id(perm.permute(v));
        }
        return permutation::Permutation(sigma);
    }

    std::vector<int> retrieve_vertices(const std::set<Level>& items) const {
        std::vector<int> used(num_variables_, 0);

        for (auto level: items) {
            graph::Vertex v = get_vertex(num_variables_ - level);
            used[v] = 1;
        }

        return used;
    }

    graph::InternalVertexId map_to_internal_vertex_id(graph::Vertex v) const {
        return mapping_vertex_[v];
    }

private:
    graph::Vertex get_vertex(graph::InternalVertexId vid) const {
        return vertex_order_[vid];
    }
};


void prepare_binary_derivative_structures_with_sro_(
    tdzdd::DdStructure<2>& dd,
    int num_sites,
    int num_types,
    const std::vector<graph::Vertex>& vertex_order,
    const std::vector<permutation::Permutation>& automorphism,
    const std::vector<permutation::Permutation>& translations,
    const std::vector<std::pair<std::vector<int>, int>>& composition_constraints,
    bool remove_superperiodic
)
{
    size_t num_variables = num_sites;
    BinaryVertexConverter converter(num_sites, vertex_order);

    // ==== translate permutations with vertex_order ====
    std::vector<permutation::Permutation> reordered_automorphism;
    reordered_automorphism.reserve(automorphism.size());
    for (const auto& perm: automorphism) {
        auto reordered_perm = converter.reorder_premutation(perm);
        reordered_automorphism.emplace_back(reordered_perm);
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
        auto reordered_perm = converter.reorder_premutation(perm);
        reordered_translation_group.emplace_back(reordered_perm);
    }

    // ======== construct DD ========
    dd = universe::Universe(num_variables);

    // spec for composition constraints
    assert(!composition_constraints.empty());
    for (const auto& variables_count: composition_constraints) {
        std::vector<int> group;
        group.reserve(variables_count.first.size());
        for (graph::Vertex v: variables_count.first) {
            group.emplace_back(converter.map_to_internal_vertex_id(v));
        }

        int k = variables_count.second; // take k elements of 1-branchs
        choice::Choice spec(num_variables, k, group, false);

        dd.zddSubset(spec);
        dd.zddReduce();
    }

    // remove symmetry duplicates
    for (const auto& perm: reordered_automorphism) {
        permutation::PermutationFrontierManager pfm(perm);
        permutation::isomorphism::IsomorphismElimination spec(pfm);
        dd.zddSubset(spec);
        dd.zddReduce();
    }

    // remove superperiodic
    if (remove_superperiodic) {
        auto identity = permutation::get_identity(num_variables);
        for (const auto& perm: reordered_translation_group) {
            if (perm == identity) {
                continue;
            }
            permutation::PermutationFrontierManager pfm(perm);
            permutation::superperiodic::SuperperiodicElimination spec(pfm);
            dd.zddSubset(spec);
            dd.zddReduce();
        }
    }

}

/// @brief enumerate DD for derivative structures with fixed short-range order (SRO)
/// @param[out] dd DD for output
/// @param[in] num_sites the number of sites in supercell
/// @param[in] num_types the kinds of species
/// @param[in] automorphism symmetry group of supercell
/// @param[in] translations (Optional) permutation group derived from
////           translations of cell. Required if `remove_superperiodic` is true.
/// @param[in] composition_constraints composition_constraints[i]
///            is a pair of sites and a desired number of label=1
/// @param[in] remove_superperiodic iff true, remove superperiodic structures
void prepare_derivative_structures_with_sro(
    tdzdd::DdStructure<2>& dd,
    int num_sites,
    int num_types,
    const std::vector<graph::Vertex>& vertex_order,
    const std::vector<permutation::Permutation>& automorphism,
    const std::vector<permutation::Permutation>& translations,
    const std::vector<std::pair<std::vector<int>, int>>& composition_constraints,
    bool remove_superperiodic
)
{
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
    // TODO: sanity check on composition_constraints
    if (remove_superperiodic && translations.empty()) {
        std::cerr << "Translational group is required for removing superperiodic structures.";
        exit(1);
    }

    tdzdd::MessageHandler::showMessages(true);
    if (num_types == 2) {
        prepare_binary_derivative_structures_with_sro_(
            dd,
            num_sites,
            num_types,
            vertex_order,
            automorphism,
            translations,
            composition_constraints,
            remove_superperiodic
        );
    } else {
        std::cerr << "TODO: multicomponent system" << std::endl;
        exit(1);
    }
    tdzdd::MessageHandler::showMessages(false);
}


/// @brief restrict a pair correlation for derivative structures
/// @param[out] dd DD for output
/// @param[in] num_sites the number of sites in supercell
/// @param[in] num_types the kinds of species
/// @param[in] graph frontier manager of cluster graph
/// @param[in] target target weight without loop offset of vertex-induced subgraph
void restrict_pair_correlation(
    tdzdd::DdStructure<2>& dd,
    int num_sites,
    int num_types,
    const std::vector<graph::Vertex>& vertex_order,
    const graph::Graph& graph,
    const graph::Weight target
)
{
    tdzdd::MessageHandler::showMessages(true);
    if (num_types == 2) {
        // need to use the same vertex-order!
        graph::VertexGraphFrontierManager vgfm(graph, vertex_order);
        graph::induced_subgraph::VertexInducedSubgraphSpec spec(vgfm, target);
        dd.zddSubset(spec);
        dd.zddReduce();
    } else {
        std::cerr << "TODO: multicomponent system" << std::endl;
        exit(1);
    }
    tdzdd::MessageHandler::showMessages(false);
}


/// @brief enumerate DD for derivative structures with fixed short-range order (SRO)
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
/// @param[in] remove_superperiodic iff true, remove superperiodic structures
void construct_derivative_structures_with_sro(
    tdzdd::DdStructure<2>& dd,
    int num_sites,
    int num_types,
    const std::vector<graph::Vertex>& vertex_order,
    const std::vector<permutation::Permutation>& automorphism,
    const std::vector<permutation::Permutation>& translations,
    const std::vector<std::pair<std::vector<int>, int>>& composition_constraints,
    const std::vector<graph::Graph>& graphs,
    const std::vector<graph::Weight>& targets,
    bool remove_superperiodic
)
{
    // ==== construct DD ====
    prepare_derivative_structures_with_sro(
        dd,
        num_sites,
        num_types,
        vertex_order,
        automorphism,
        translations,
        composition_constraints,
        remove_superperiodic
    );

    // fix SRO
    // TODO: better vertex_order choice
    size_t num_graphs = graphs.size();
    assert(targets.size() == num_graphs);
    for (size_t i = 0; i < num_graphs; ++i) {
        // need to use the same vertex-order!
        restrict_pair_correlation(
            dd,
            num_sites,
            num_types,
            vertex_order,
            graphs[i],
            targets[i]
        );
    }
}


std::vector<int> convert_to_binary_labeling_with_graph(
    const tdzdd::DdStructure<2>::const_iterator& itr,
    const BinaryVertexConverter& converter
)
{
    // vertex order in DD is diffenrent from the original variable order in the graph
    std::vector<int> labeling = converter.retrieve_vertices(*itr);
    return labeling;
}

} // namespace derivative_structure
} // namespace pyzdd

#endif // PYZDD_SHORT_RANGE_ORDRE_H
