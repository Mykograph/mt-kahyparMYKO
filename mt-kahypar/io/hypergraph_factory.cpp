/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2023 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#include "hypergraph_factory.h"

#include <fstream>
#include <limits>

#include "mt-kahypar/macros.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/datastructures/fixed_vertex_support.h"
#include "mt-kahypar/partition/conversion.h"
#include "mt-kahypar/utils/exception.h"

namespace mt_kahypar {
namespace io {

// OriginalEdgeSnapshot is declared in hypergraph_factory.h

namespace {

template<typename Hypergraph>
mt_kahypar_hypergraph_t constructHypergraph(const HypernodeID& num_hypernodes,
                                            const HyperedgeID& num_hyperedges,
                                            const HyperedgeVector& hyperedges,
                                            const HyperedgeWeight* hyperedge_weight,
                                            const HypernodeWeight* hypernode_weight,
                                            const HypernodeID num_removed_single_pin_hes,
                                            const bool stable_construction) {
  Hypergraph* hypergraph = new Hypergraph();
  *hypergraph = Hypergraph::Factory::construct(num_hypernodes, num_hyperedges, hyperedges,
    hyperedge_weight, hypernode_weight, stable_construction);
  hypergraph->setNumRemovedHyperedges(num_removed_single_pin_hes);
  return mt_kahypar_hypergraph_t {
    reinterpret_cast<mt_kahypar_hypergraph_s*>(hypergraph), Hypergraph::TYPE };
}

mt_kahypar_hypergraph_t constructHypergraph(const mt_kahypar_hypergraph_type_t& type,
                                            const HypernodeID& num_hypernodes,
                                            const HyperedgeID& num_hyperedges,
                                            const HyperedgeVector& hyperedges,
                                            vec<HyperedgeWeight>& hyperedge_weight,
                                            vec<HypernodeWeight>& hypernode_weight,
                                            const HypernodeID num_removed_single_pin_hes,
                                            const bool stable_construction) {
  switch ( type ) {
    case STATIC_HYPERGRAPH:
      return constructHypergraph<ds::StaticHypergraph>(
        num_hypernodes, num_hyperedges, hyperedges,
        hyperedge_weight.data(), hypernode_weight.data(),
        num_removed_single_pin_hes, stable_construction);
    case STATIC_GRAPH:
      ENABLE_GRAPHS(
        return constructHypergraph<ds::StaticGraph>(
          num_hypernodes, num_hyperedges, hyperedges,
          hyperedge_weight.data(), hypernode_weight.data(),
          num_removed_single_pin_hes, stable_construction);
      )
    case DYNAMIC_HYPERGRAPH:
      ENABLE_HIGHEST_QUALITY(
        return constructHypergraph<ds::DynamicHypergraph>(
          num_hypernodes, num_hyperedges, hyperedges,
          hyperedge_weight.data(), hypernode_weight.data(),
          num_removed_single_pin_hes, stable_construction);
      )
    case DYNAMIC_GRAPH:
      ENABLE_HIGHEST_QUALITY_FOR_GRAPHS(
        return constructHypergraph<ds::DynamicGraph>(
          num_hypernodes, num_hyperedges, hyperedges,
          hyperedge_weight.data(), hypernode_weight.data(),
          num_removed_single_pin_hes, stable_construction);
      )
    case NULLPTR_HYPERGRAPH:
      return mt_kahypar_hypergraph_t { nullptr, NULLPTR_HYPERGRAPH };
  }
  return mt_kahypar_hypergraph_t { nullptr, NULLPTR_HYPERGRAPH };
}

template<typename Hypergraph>
mt_kahypar_hypergraph_t constructGraph(const HypernodeID& num_nodes,
                                       const HyperedgeID& num_edges,
                                       const EdgeVector& edges,
                                       const HyperedgeWeight* edge_weight,
                                       const HypernodeWeight* node_weight,
                                       const bool stable_construction) {
  Hypergraph* graph = new Hypergraph();
  *graph = Hypergraph::Factory::construct_from_graph_edges(num_nodes, num_edges, edges,
    edge_weight, node_weight, stable_construction);
  return mt_kahypar_hypergraph_t {
    reinterpret_cast<mt_kahypar_hypergraph_s*>(graph), Hypergraph::TYPE };
}

mt_kahypar_hypergraph_t constructGraph(const mt_kahypar_hypergraph_type_t& type,
                                       const HypernodeID& num_nodes,
                                       const HyperedgeID& num_edges,
                                       const EdgeVector& edges,
                                       vec<HyperedgeWeight>& edge_weight,
                                       vec<HypernodeWeight>& node_weight,
                                       const bool stable_construction) {
  switch ( type ) {
    case STATIC_GRAPH:
      ENABLE_GRAPHS(
        return constructGraph<ds::StaticGraph>(
          num_nodes, num_edges, edges,
          edge_weight.data(), node_weight.data(), stable_construction);
      )
    case DYNAMIC_GRAPH:
      ENABLE_HIGHEST_QUALITY_FOR_GRAPHS(
        return constructGraph<ds::DynamicGraph>(
          num_nodes, num_edges, edges,
          edge_weight.data(), node_weight.data(), stable_construction);
      )
    case STATIC_HYPERGRAPH:
    case DYNAMIC_HYPERGRAPH:
    case NULLPTR_HYPERGRAPH:
      return mt_kahypar_hypergraph_t { nullptr, NULLPTR_HYPERGRAPH };
  }
  return mt_kahypar_hypergraph_t { nullptr, NULLPTR_HYPERGRAPH };
}

struct ConstraintPairHash {
  size_t operator()(const std::pair<HypernodeID, HypernodeID>& p) const noexcept {
    return (static_cast<size_t>(p.first) << 32) ^ static_cast<size_t>(p.second);
  }
};

// Overload for HyperedgeVector (hypergraph case) — edges accessed via edge[0], edge[1]
//
// `normal_edge_weight_bonus` (x) is added to every edge weight that is NOT an
// anti-constraint edge (i.e. every "normal" edge), whether it was already
// present in the input or newly appended by a later constraint pass call.
void applyConstraintPairs(const std::string& constraint_filename,
                          const HyperedgeWeight constraint_weight,
                          const HypernodeID num_nodes,
                          HyperedgeID& num_edges,
                          HyperedgeVector& edges,
                          vec<HyperedgeWeight>& edge_weight,
                          const HyperedgeWeight normal_edge_weight_bonus = 0,
                          vec<std::pair<HypernodeID, HypernodeID>>* out_constraints = nullptr) {
  if (constraint_filename.empty()) {
    // No constraints to apply, but the normal-edge bonus still applies to all edges.
    if (normal_edge_weight_bonus != 0) {
      for (HyperedgeID he = 0; he < static_cast<HyperedgeID>(edge_weight.size()); ++he) {
        edge_weight[he] += normal_edge_weight_bonus;
      }
    }
    return;
  }

  std::ifstream file(constraint_filename);
  if (!file) throw InvalidInputException("File not found: " + constraint_filename);

  vec<std::pair<HypernodeID, HypernodeID>> constraints;
  int64_t u = 0, v = 0;
  while (file >> u >> v) {
    if (u < 0 || v < 0)
      throw InvalidInputException("Constraint file contains a negative node id: " + constraint_filename);
    if (u >= num_nodes || v >= num_nodes)
      throw InvalidInputException("Constraint file contains an invalid node id: " + constraint_filename);
    if (u == v)
      throw InvalidInputException("Constraint file contains a self-constraint for node " + std::to_string(u));
    HypernodeID hu = static_cast<HypernodeID>(u);
    HypernodeID hv = static_cast<HypernodeID>(v);
    if (hu > hv) std::swap(hu, hv);
    constraints.emplace_back(hu, hv);
  }

  if (constraints.empty()) {
    if (normal_edge_weight_bonus != 0) {
      for (HyperedgeID he = 0; he < static_cast<HyperedgeID>(edge_weight.size()); ++he) {
        edge_weight[he] += normal_edge_weight_bonus;
      }
    }
    return;
  }

  std::unordered_map<std::pair<HypernodeID, HypernodeID>, vec<HyperedgeID>, ConstraintPairHash> edge_map;
  edge_map.reserve(edges.size());
  for (HyperedgeID he = 0; he < static_cast<HyperedgeID>(edges.size()); ++he) {
    const auto& edge = edges[he];
    if (edge.size() != 2) continue;
    HypernodeID eu = edge[0], ev = edge[1];
    if (eu > ev) std::swap(eu, ev);
    edge_map[{eu, ev}].push_back(he);
  }

  // Track which edge indices are anti-constraint edges so the normal-edge
  // bonus below can skip them.
  vec<bool> is_constraint_edge(edges.size(), false);

  for (const auto& constraint : constraints) {
    auto it = edge_map.find(constraint);
    if (it != edge_map.end()) {
      for (const HyperedgeID he : it->second) {
        edge_weight[he] = constraint_weight;
        is_constraint_edge[he] = true;
      }
    } else {
      edges.emplace_back();
      auto& new_edge = edges.back();
      new_edge.reserve(2);
      new_edge.push_back(constraint.first);
      new_edge.push_back(constraint.second);
      edge_weight.push_back(constraint_weight);
      is_constraint_edge.push_back(true);
      edge_map[constraint].push_back(static_cast<HyperedgeID>(edges.size() - 1));
      ++num_edges;
    }
  }

  // Add the preset bonus x to every edge that is NOT an anti-constraint edge.
  if (normal_edge_weight_bonus != 0) {
    for (HyperedgeID he = 0; he < static_cast<HyperedgeID>(edge_weight.size()); ++he) {
      if (!is_constraint_edge[he]) {
        edge_weight[he] += normal_edge_weight_bonus;
      }
    }
  }

  if (out_constraints) *out_constraints = constraints;
}

// Overload for EdgeVector (graph case) — edges accessed via .first, .second
//
// `normal_edge_weight_bonus` (x) is added to every edge weight that is NOT an
// anti-constraint edge (i.e. every "normal" edge), whether it was already
// present in the input or newly appended by a later constraint pass call.
void applyConstraintPairs(const std::string& constraint_filename,
                          const HyperedgeWeight constraint_weight,
                          const HypernodeID num_nodes,
                          HyperedgeID& num_edges,
                          EdgeVector& edges,
                          vec<HyperedgeWeight>& edge_weight,
                          const HyperedgeWeight normal_edge_weight_bonus = 0,
                          vec<std::pair<HypernodeID, HypernodeID>>* out_constraints = nullptr) {
  if (constraint_filename.empty()) {
    if (normal_edge_weight_bonus != 0) {
      for (HyperedgeID he = 0; he < static_cast<HyperedgeID>(edge_weight.size()); ++he) {
        edge_weight[he] += normal_edge_weight_bonus;
      }
    }
    return;
  }

  std::ifstream file(constraint_filename);
  if (!file) throw InvalidInputException("File not found: " + constraint_filename);

  vec<std::pair<HypernodeID, HypernodeID>> constraints;
  int64_t u = 0, v = 0;
  while (file >> u >> v) {
    if (u < 0 || v < 0)
      throw InvalidInputException("Constraint file contains a negative node id: " + constraint_filename);
    if (u >= num_nodes || v >= num_nodes)
      throw InvalidInputException("Constraint file contains an invalid node id: " + constraint_filename);
    if (u == v)
      throw InvalidInputException("Constraint file contains a self-constraint for node " + std::to_string(u));
    HypernodeID hu = static_cast<HypernodeID>(u);
    HypernodeID hv = static_cast<HypernodeID>(v);
    if (hu > hv) std::swap(hu, hv);
    constraints.emplace_back(hu, hv);
  }

  if (constraints.empty()) {
    if (normal_edge_weight_bonus != 0) {
      for (HyperedgeID he = 0; he < static_cast<HyperedgeID>(edge_weight.size()); ++he) {
        edge_weight[he] += normal_edge_weight_bonus;
      }
    }
    return;
  }

  std::unordered_map<std::pair<HypernodeID, HypernodeID>, vec<HyperedgeID>, ConstraintPairHash> edge_map;
  edge_map.reserve(edges.size());
  for (HyperedgeID he = 0; he < static_cast<HyperedgeID>(edges.size()); ++he) {
    HypernodeID eu = edges[he].first, ev = edges[he].second;
    if (eu > ev) std::swap(eu, ev);
    edge_map[{eu, ev}].push_back(he);
  }

  vec<bool> is_constraint_edge(edges.size(), false);

  for (const auto& constraint : constraints) {
    auto it = edge_map.find(constraint);
    if (it != edge_map.end()) {
      for (const HyperedgeID he : it->second) {
        edge_weight[he] = constraint_weight;
        is_constraint_edge[he] = true;
      }
    } else {
      edges.emplace_back(constraint.first, constraint.second);
      edge_weight.push_back(constraint_weight);
      is_constraint_edge.push_back(true);
      edge_map[constraint].push_back(static_cast<HyperedgeID>(edges.size() - 1));
      ++num_edges;
    }
  }

  if (normal_edge_weight_bonus != 0) {
    for (HyperedgeID he = 0; he < static_cast<HyperedgeID>(edge_weight.size()); ++he) {
      if (!is_constraint_edge[he]) {
        edge_weight[he] += normal_edge_weight_bonus;
      }
    }
  }

  if (out_constraints) *out_constraints = constraints;
}

// Helper to attach the constraint graph to a typed hypergraph after construction
template<typename Hypergraph>
void attachConstraintGraph(Hypergraph& hg,
                           const vec<std::pair<HypernodeID, HypernodeID>>& constraints) {
  hg.fixedVertexSupport().setHypergraph(&hg);
  hg.fixedVertexSupport().setNegativeConstraints(constraints);
}

mt_kahypar_hypergraph_t readHMetisFile(const std::string& filename,
                                       const mt_kahypar_hypergraph_type_t& type,
                                       const bool stable_construction,
                                       const bool remove_single_pin_hes,
                                       const bool print_warnings,
                                       const std::string& constraint_filename,
                                       const HyperedgeWeight constraint_weight,
                                       const HyperedgeWeight normal_edge_weight_bonus = 0,
                                       OriginalEdgeSnapshot* out_snapshot = nullptr) {
  HyperedgeID num_hyperedges = 0;
  HypernodeID num_hypernodes = 0;
  HyperedgeID num_removed_single_pin_hyperedges = 0;
  HyperedgeVector hyperedges;
  vec<HyperedgeWeight> hyperedges_weight;
  vec<HypernodeWeight> hypernodes_weight;
  readHypergraphFile(filename, num_hyperedges, num_hypernodes,
                     num_removed_single_pin_hyperedges, hyperedges,
                     hyperedges_weight, hypernodes_weight, remove_single_pin_hes, print_warnings);

  // Snapshot BEFORE the constraint pass mutates weights or appends edges.
  if (out_snapshot) {
    out_snapshot->is_graph = false;
    out_snapshot->valid = true;
    out_snapshot->original_hyperedges = hyperedges;          // deep copy
    out_snapshot->original_edge_weights = hyperedges_weight; // deep copy
  }

  vec<std::pair<HypernodeID, HypernodeID>> constraints;
  applyConstraintPairs(constraint_filename, constraint_weight, num_hypernodes,
                       num_hyperedges, hyperedges, hyperedges_weight,
                       normal_edge_weight_bonus, &constraints);

  mt_kahypar_hypergraph_t hg = constructHypergraph(type, num_hypernodes, num_hyperedges,
                                                    hyperedges, hyperedges_weight,
                                                    hypernodes_weight,
                                                    num_removed_single_pin_hyperedges,
                                                    stable_construction);
  if (!constraints.empty()) {
    switch (hg.type) {
      case STATIC_HYPERGRAPH:
        attachConstraintGraph(utils::cast<ds::StaticHypergraph>(hg), constraints); break;
      ENABLE_HIGHEST_QUALITY(case DYNAMIC_HYPERGRAPH:
        attachConstraintGraph(utils::cast<ds::DynamicHypergraph>(hg), constraints); break;)
      default: break;
    }
  }
  return hg;
}

mt_kahypar_hypergraph_t readMetisFile(const std::string& filename,
                                      const mt_kahypar_hypergraph_type_t& type,
                                      const bool stable_construction,
                                      const bool,
                                      const std::string& constraint_filename,
                                      const HyperedgeWeight constraint_weight,
                                      const HyperedgeWeight normal_edge_weight_bonus = 0,
                                      OriginalEdgeSnapshot* out_snapshot = nullptr) {
  HyperedgeID num_edges = 0;
  HypernodeID num_vertices = 0;
  vec<HyperedgeWeight> edges_weight;
  vec<HypernodeWeight> nodes_weight;

  if (type == STATIC_GRAPH || type == DYNAMIC_GRAPH) {
    EdgeVector edges;
    vec<std::pair<HypernodeID, HypernodeID>> constraints;
    readGraphFile(filename, num_edges, num_vertices, edges, edges_weight, nodes_weight);

    // Snapshot BEFORE the constraint pass mutates weights or appends edges.
    if (out_snapshot) {
      out_snapshot->is_graph = true;
      out_snapshot->valid = true;
      out_snapshot->original_edges = edges;               // deep copy
      out_snapshot->original_edge_weights = edges_weight; // deep copy
    }

    applyConstraintPairs(constraint_filename, constraint_weight, num_vertices,
                         num_edges, edges, edges_weight,
                         normal_edge_weight_bonus, &constraints);
    mt_kahypar_hypergraph_t hg = constructGraph(type, num_vertices, num_edges, edges,
                                                edges_weight, nodes_weight, stable_construction);
    if (!constraints.empty()) {
      switch (hg.type) {
        ENABLE_GRAPHS(case STATIC_GRAPH:
          attachConstraintGraph(utils::cast<ds::StaticGraph>(hg), constraints); break;)
        ENABLE_HIGHEST_QUALITY_FOR_GRAPHS(case DYNAMIC_GRAPH:
          attachConstraintGraph(utils::cast<ds::DynamicGraph>(hg), constraints); break;)
        default: break;
      }
    }
    return hg;
  } else {
    HyperedgeVector edges;
    readGraphFile(filename, num_edges, num_vertices, edges, edges_weight, nodes_weight);

    // Snapshot BEFORE the constraint pass mutates weights or appends edges.
    if (out_snapshot) {
      out_snapshot->is_graph = false;
      out_snapshot->valid = true;
      out_snapshot->original_hyperedges = edges;          // deep copy
      out_snapshot->original_edge_weights = edges_weight; // deep copy
    }

    applyConstraintPairs(constraint_filename, constraint_weight, num_vertices,
                         num_edges, edges, edges_weight, normal_edge_weight_bonus);
    return constructHypergraph(type, num_vertices, num_edges, edges,
                               edges_weight, nodes_weight, 0, stable_construction);
  }
}

} // namespace

mt_kahypar_hypergraph_t readInputFile(const std::string& filename,
                                      const PresetType& preset,
                                      const InstanceType& instance,
                                      const FileFormat& format,
                                      const bool stable_construction,
                                      const bool remove_single_pin_hes,
                                      const bool print_warnings,
                                      const std::string& constraint_filename,
                                      const HyperedgeWeight constraint_weight,
                                      const HyperedgeWeight normal_edge_weight_bonus,
                                      OriginalEdgeSnapshot* out_snapshot) {
  mt_kahypar_hypergraph_type_t type = to_hypergraph_c_type(preset, instance);
  switch ( format ) {
    case FileFormat::hMetis: return readHMetisFile(
      filename, type, stable_construction, remove_single_pin_hes, print_warnings,
      constraint_filename, constraint_weight, normal_edge_weight_bonus, out_snapshot);
    case FileFormat::Metis: return readMetisFile(
      filename, type, stable_construction, print_warnings,
      constraint_filename, constraint_weight, normal_edge_weight_bonus, out_snapshot);
  }
  return mt_kahypar_hypergraph_t { nullptr, NULLPTR_HYPERGRAPH };
}

template<typename Hypergraph>
Hypergraph readInputFile(const std::string& filename,
                         const FileFormat& format,
                         const bool stable_construction,
                         const bool remove_single_pin_hes,
                         const bool print_warnings,
                         const std::string& constraint_filename,
                         const HyperedgeWeight constraint_weight,
                         const HyperedgeWeight normal_edge_weight_bonus,
                         OriginalEdgeSnapshot* out_snapshot) {
  mt_kahypar_hypergraph_t hypergraph { nullptr, NULLPTR_HYPERGRAPH };
  switch ( format ) {
    case FileFormat::hMetis: hypergraph = readHMetisFile(
      filename, Hypergraph::TYPE, stable_construction, remove_single_pin_hes, print_warnings,
      constraint_filename, constraint_weight, normal_edge_weight_bonus, out_snapshot);
      break;
    case FileFormat::Metis: hypergraph = readMetisFile(
      filename, Hypergraph::TYPE, stable_construction, print_warnings,
      constraint_filename, constraint_weight, normal_edge_weight_bonus, out_snapshot);
  }
  return std::move(utils::cast<Hypergraph>(hypergraph));
}

// ----------------------------------------------------------------------------
// Cut counting against the ORIGINAL (un-tampered) edge set captured in the
// snapshot. `part_id` must map HypernodeID -> PartitionID for the final
// partition (the same node IDs are used before and after the constraint
// pass, since constraints only add edges / change weights, never nodes).
// ----------------------------------------------------------------------------

HyperedgeWeight countOriginalCutHyperedges(const HyperedgeVector& original_hyperedges,
                                           const vec<HyperedgeWeight>& original_weights,
                                           const vec<PartitionID>& part_id) {
  HyperedgeWeight cut = 0;
  for (size_t he = 0; he < original_hyperedges.size(); ++he) {
    const auto& pins = original_hyperedges[he];
    if (pins.empty()) continue;
    const PartitionID first_block = part_id[pins[0]];
    for (const HypernodeID pin : pins) {
      if (part_id[pin] != first_block) {
        cut += original_weights[he];
        break;
      }
    }
  }
  return cut;
}

HyperedgeWeight countOriginalCutEdges(const EdgeVector& original_edges,
                                      const vec<HyperedgeWeight>& original_weights,
                                      const vec<PartitionID>& part_id) {
  HyperedgeWeight cut = 0;
  for (size_t e = 0; e < original_edges.size(); ++e) {
    if (part_id[original_edges[e].first] != part_id[original_edges[e].second]) {
      cut += original_weights[e];
    }
  }
  return cut;
}

HyperedgeWeight countOriginalCut(const OriginalEdgeSnapshot& snapshot,
                                 const vec<PartitionID>& part_id) {
  if (!snapshot.valid) {
    throw InvalidInputException("Attempted to count cut on an empty/uninitialized OriginalEdgeSnapshot");
  }
  if (snapshot.is_graph) {
    return countOriginalCutEdges(snapshot.original_edges, snapshot.original_edge_weights, part_id);
  } else {
    return countOriginalCutHyperedges(snapshot.original_hyperedges, snapshot.original_edge_weights, part_id);
  }
}

namespace {

HypernodeID numberOfNodes(mt_kahypar_hypergraph_t hypergraph) {
  switch ( hypergraph.type ) {
    case STATIC_HYPERGRAPH: return utils::cast<ds::StaticHypergraph>(hypergraph).initialNumNodes();
    ENABLE_GRAPHS(case STATIC_GRAPH: return utils::cast<ds::StaticGraph>(hypergraph).initialNumNodes();)
    ENABLE_HIGHEST_QUALITY(case DYNAMIC_HYPERGRAPH: return utils::cast<ds::DynamicHypergraph>(hypergraph).initialNumNodes();)
    ENABLE_HIGHEST_QUALITY_FOR_GRAPHS(case DYNAMIC_GRAPH: return utils::cast<ds::DynamicGraph>(hypergraph).initialNumNodes();)
    case NULLPTR_HYPERGRAPH: return 0;
    default: return 0;
  }
}

template<typename Hypergraph>
void addFixedVertices(Hypergraph& hypergraph,
                      const mt_kahypar_partition_id_t* fixed_vertices,
                      const PartitionID k) {
  ds::FixedVertexSupport<Hypergraph> fixed_vertex_support(
    hypergraph.initialNumNodes(), k);
  fixed_vertex_support.setHypergraph(&hypergraph);
  hypergraph.doParallelForAllNodes([&](const HypernodeID& hn) {
    if ( fixed_vertices[hn] != -1 ) {
      if ( fixed_vertices[hn] < 0 || fixed_vertices[hn] >= k ) {
        throw InvalidInputException(
          "Try to partition hypergraph into " + STR(k) + " blocks, but node " +
           STR(hn) + " is fixed to block " + STR(fixed_vertices[hn]));
      }
      fixed_vertex_support.fixToBlock(hn, fixed_vertices[hn]);
    }
  });
  hypergraph.addFixedVertexSupport(std::move(fixed_vertex_support));
}

template<typename Hypergraph>
void removeFixedVertices(Hypergraph& hypergraph) {
  ds::FixedVertexSupport<Hypergraph> fixed_vertex_support;
  hypergraph.addFixedVertexSupport(std::move(fixed_vertex_support));
}

} // namespace

void addFixedVertices(mt_kahypar_hypergraph_t hypergraph,
                      const mt_kahypar_partition_id_t* fixed_vertices,
                      const PartitionID k) {
  switch ( hypergraph.type ) {
    case STATIC_HYPERGRAPH:
      addFixedVertices(utils::cast<ds::StaticHypergraph>(hypergraph), fixed_vertices, k); break;
    ENABLE_GRAPHS(case STATIC_GRAPH:
      addFixedVertices(utils::cast<ds::StaticGraph>(hypergraph), fixed_vertices, k); break;)
    ENABLE_HIGHEST_QUALITY(case DYNAMIC_HYPERGRAPH:
      addFixedVertices(utils::cast<ds::DynamicHypergraph>(hypergraph), fixed_vertices, k); break;)
    ENABLE_HIGHEST_QUALITY_FOR_GRAPHS(case DYNAMIC_GRAPH:
      addFixedVertices(utils::cast<ds::DynamicGraph>(hypergraph), fixed_vertices, k); break;)
    case NULLPTR_HYPERGRAPH:
    default: break;
  }
}

void addFixedVerticesFromFile(mt_kahypar_hypergraph_t hypergraph,
                              const std::string& filename,
                              const PartitionID k) {
  std::vector<PartitionID> fixed_vertices;
  io::readPartitionFile(filename, numberOfNodes(hypergraph), fixed_vertices);
  addFixedVertices(hypergraph, fixed_vertices.data(), k);
}

void removeFixedVertices(mt_kahypar_hypergraph_t hypergraph) {
  switch ( hypergraph.type ) {
    case STATIC_HYPERGRAPH:
      removeFixedVertices(utils::cast<ds::StaticHypergraph>(hypergraph)); break;
    ENABLE_GRAPHS(case STATIC_GRAPH:
      removeFixedVertices(utils::cast<ds::StaticGraph>(hypergraph)); break;)
    ENABLE_HIGHEST_QUALITY(case DYNAMIC_HYPERGRAPH:
      removeFixedVertices(utils::cast<ds::DynamicHypergraph>(hypergraph)); break;)
    ENABLE_HIGHEST_QUALITY_FOR_GRAPHS(case DYNAMIC_GRAPH:
      removeFixedVertices(utils::cast<ds::DynamicGraph>(hypergraph)); break;)
    case NULLPTR_HYPERGRAPH:
    default: break;
  }
}

template ds::StaticHypergraph readInputFile(const std::string& filename,
                                            const FileFormat& format,
                                            const bool stable_construction,
                                            const bool remove_single_pin_hes,
                                            const bool logging,
                                            const std::string& constraint_filename,
                                            const HyperedgeWeight constraint_weight,
                                            const HyperedgeWeight normal_edge_weight_bonus,
                                            OriginalEdgeSnapshot* out_snapshot);

ENABLE_GRAPHS(template ds::StaticGraph readInputFile(const std::string& filename,
                                                     const FileFormat& format,
                                                     const bool stable_construction,
                                                     const bool remove_single_pin_hes,
                                                     const bool logging,
                                                     const std::string& constraint_filename,
                                                     const HyperedgeWeight constraint_weight,
                                                     const HyperedgeWeight normal_edge_weight_bonus,
                                                     OriginalEdgeSnapshot* out_snapshot);)

ENABLE_HIGHEST_QUALITY(template ds::DynamicHypergraph readInputFile(const std::string& filename,
                                                                    const FileFormat& format,
                                                                    const bool stable_construction,
                                                                    const bool remove_single_pin_hes,
                                                                    const bool logging,
                                                                    const std::string& constraint_filename,
                                                                    const HyperedgeWeight constraint_weight,
                                                                    const HyperedgeWeight normal_edge_weight_bonus,
                                                                    OriginalEdgeSnapshot* out_snapshot);)

ENABLE_HIGHEST_QUALITY_FOR_GRAPHS(template ds::DynamicGraph readInputFile(const std::string& filename,
                                                                          const FileFormat& format,
                                                                          const bool stable_construction,
                                                                          const bool remove_single_pin_hes,
                                                                          const bool logging,
                                                                          const std::string& constraint_filename,
                                                                          const HyperedgeWeight constraint_weight,
                                                                          const HyperedgeWeight normal_edge_weight_bonus,
                                                                          OriginalEdgeSnapshot* out_snapshot);)

#ifndef KAHYPAR_ENABLE_GRAPH_PARTITIONING_FEATURES
template ds::StaticGraph readInputFile(const std::string& filename,
                                       const FileFormat& format,
                                       const bool stable_construction,
                                       const bool remove_single_pin_hes,
                                       const bool logging,
                                       const std::string& constraint_filename,
                                       const HyperedgeWeight constraint_weight,
                                       const HyperedgeWeight normal_edge_weight_bonus,
                                       OriginalEdgeSnapshot* out_snapshot);
#endif

}  // namespace io
}  // namespace mt_kahypar