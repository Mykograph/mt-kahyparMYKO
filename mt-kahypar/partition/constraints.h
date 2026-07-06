#pragma once

#include <fstream>
#include <limits>
#include <memory>
#include <vector>

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/factories.h"
#include "mt-kahypar/partition/refinement/gains/gain_cache_ptr.h"
#include "mt-kahypar/datastructures/priority_queue.h"
#include "mt-kahypar/utils/cast.h"

namespace mt_kahypar::constraints {

// ####################### Helpers #######################

// Reads constraint pairs directly from file — no dependency on io:: functions
inline vec<std::pair<HypernodeID, HypernodeID>> readConstraintFile(const std::string& filename) {
  vec<std::pair<HypernodeID, HypernodeID>> constraints;
  if (filename.empty()) return constraints;
  std::ifstream file(filename);
  if (!file) return constraints;
  int64_t u = 0, v = 0;
  while (file >> u >> v) {
    if (u >= 0 && v >= 0 && u != v) {
      constraints.emplace_back(
        static_cast<HypernodeID>(u),
        static_cast<HypernodeID>(v));
    }
  }
  return constraints;
}

// ####################### Verify #######################

template<typename PartitionedHypergraph>
bool verifyConstraints(const PartitionedHypergraph& partitioned_hg, const Context& context) {
  const auto constraints = readConstraintFile(context.partition.constraint_filename);
  for (const auto& [u, v] : constraints) {
    if (partitioned_hg.partID(u) == partitioned_hg.partID(v)) {
      return false;
    }
  }
  return true;
}

template<typename PartitionedHypergraph>
size_t countConstraintViolations(const PartitionedHypergraph& partitioned_hg, const Context& context) {
  const auto constraints = readConstraintFile(context.partition.constraint_filename);
  size_t cnt = 0;
  for (const auto& [u, v] : constraints) {
    if (partitioned_hg.partID(u) == partitioned_hg.partID(v)) ++cnt;
  }
  return cnt;
}

// ####################### PQ setup #######################

using Key = std::tuple<HypernodeID, HypernodeID, HypernodeID>;

template<typename Comparator = std::less<Key>, uint32_t arity = 4>
using PQ = ds::Heap<Key, HypernodeID, Comparator, arity>;

// ####################### Core PP functions #######################

// Counts how many constraint-graph neighbors of a constraint node are in the same block
template<typename PartitionedHypergraph>
HypernodeID incidentNodesInSamePart(const PartitionedHypergraph& partitioned_hg,
                                    const HypernodeID constraint_node_id) {
  const ds::DynamicGraph& cg = partitioned_hg.fixedVertexSupport().getConstraintGraph();
  HypernodeID hg_node_id = static_cast<HypernodeID>(cg.nodeWeight(constraint_node_id));

  if (partitioned_hg.partID(hg_node_id) == kInvalidPartition) return 0;
  PartitionID part = partitioned_hg.partID(hg_node_id);

  HypernodeID num_nodes = 0;
  for (const auto& edge_id : cg.incidentEdges(constraint_node_id)) {
    HypernodeID neighbor = cg.edge(edge_id).target;
    HypernodeID neighbor_hg_id = static_cast<HypernodeID>(cg.nodeWeight(neighbor));
    if (partitioned_hg.partID(neighbor_hg_id) == part) num_nodes++;
  }
  return num_nodes;
}

// Picks the valid partition with the highest gain; falls back to current if all invalid
template<typename PartitionedHypergraph>
PartitionID getBestPartition(const HypernodeID node_id,
                             const PartitionID current_partition,
                             const vec<bool>& is_partition_invalid,
                             gain_cache_t& gain_cache) {
  PartitionID best_partition = current_partition;
  const PartitionID num_partitions = static_cast<PartitionID>(is_partition_invalid.size());
  HyperedgeWeight max_gain = std::numeric_limits<HyperedgeWeight>::min();

  for (PartitionID p = 0; p < num_partitions; p++) {
    if (is_partition_invalid[p]) continue;
    HyperedgeWeight gain = GainCachePtr::applyWithConcreteGainCacheForHG<PartitionedHypergraph>(
      [&](const auto& cache) {
        return cache.gain(node_id, current_partition, p);
      },
      gain_cache
    );
    if (gain > max_gain) {
      max_gain = gain;
      best_partition = p;
    }
  }
  return best_partition;
}

// Processes nodes in descending order of violation count, moving them to fix constraints
template<typename PartitionedHypergraph>
void descendingConstraintDegree(PartitionedHypergraph& partitioned_hg,
                                gain_cache_t& gain_cache) {
  GainCachePtr::applyWithConcreteGainCacheForHG<PartitionedHypergraph>(
    [&](auto& cache) {
      cache.initializeGainCache(partitioned_hg);
    },
    gain_cache
  );

  auto delta_func = [&partitioned_hg, &gain_cache](const SynchronizedEdgeUpdate& sync_update) {
    GainCachePtr::applyWithConcreteGainCacheForHG<PartitionedHypergraph>(
      [&](auto& cache) {
        cache.deltaGainUpdate(partitioned_hg, sync_update);
      },
      gain_cache
    );
  };

  const ds::DynamicGraph& cg = 
  const_cast<const PartitionedHypergraph&>(partitioned_hg).fixedVertexSupport().getConstraintGraph();
 std::vector<PosT> positions(cg.initialNumNodes(), invalid_position);
  PQ<std::less<Key>, 4> heap(positions.data(), positions.size());

  // Insert all constraint graph nodes, keyed by (violations, degree, node_id)
  for (auto node : cg.nodes()) {
    heap.insert(node, {
      incidentNodesInSamePart(partitioned_hg, node),
      cg.nodeDegree(node),
      node
    });
  }

  while (!heap.empty()) {
    HypernodeID cnode = heap.top();
    heap.deleteTop();

    // Recover the original hypergraph node ID from the constraint graph node weight
    HypernodeID hg_node_id = static_cast<HypernodeID>(cg.nodeWeight(cnode));
    PartitionID partition_id = partitioned_hg.partID(hg_node_id);

    // Mark all partitions that neighbors currently occupy as invalid
    vec<bool> invalid_partitions(partitioned_hg.k(), false);
    for (const auto& edge_id : cg.incidentEdges(cnode)) {
      HypernodeID neighbor = cg.edge(edge_id).target;
      HypernodeID neighbor_hg_id = static_cast<HypernodeID>(cg.nodeWeight(neighbor));
      invalid_partitions[partitioned_hg.partID(neighbor_hg_id)] = true;
    }

    // Only move if currently violating
    if (!invalid_partitions[partition_id]) continue;

    PartitionID new_partition = getBestPartition<PartitionedHypergraph>(
      hg_node_id, partition_id, invalid_partitions, gain_cache);

    if (new_partition != partition_id) {
      partitioned_hg.changeNodePart(hg_node_id, partition_id, new_partition, delta_func);

      // Update keys of all constraint-graph neighbors so their violation
      // counts reflect the recent move.
      for (const auto& edge_id : cg.incidentEdges(cnode)) {
        HypernodeID neighbor = cg.edge(edge_id).target;
        Key newKey = {
          incidentNodesInSamePart(partitioned_hg, neighbor),
          cg.nodeDegree(neighbor),
          neighbor
        };
        if (heap.contains(neighbor)) {
          heap.adjustKey(neighbor, newKey);
        } else if (static_cast<size_t>(neighbor) < positions.size()) {
          heap.insert(neighbor, newKey);
        }
      }

      // Re-evaluate the moved constraint node and reinsert if still violating
      Key curKey = {
        incidentNodesInSamePart(partitioned_hg, cnode),
        cg.nodeDegree(cnode),
        cnode
      };
      if (std::get<0>(curKey) > 0) {
        if (heap.contains(cnode)) {
          heap.adjustKey(cnode, curKey);
        } else if (static_cast<size_t>(cnode) < positions.size()) {
          heap.insert(cnode, curKey);
        }
      }
    }
  }
}

// ####################### Entry point #######################

template<typename PartitionedHypergraph>
void postprocessNegativeConstraints(PartitionedHypergraph& partitioned_hg,
                                    const Context& context) {
  gain_cache_t gain_cache = GainCachePtr::constructGainCache(context);

  // Repeat adaptive PQ passes until no violations remain or no improvement
  const int max_iters = 5;
  size_t prev_violations = countConstraintViolations(partitioned_hg, context);
  for (int iter = 0; iter < max_iters; ++iter) {
    descendingConstraintDegree(partitioned_hg, gain_cache);
    size_t cur_violations = countConstraintViolations(partitioned_hg, context);
    LOG << "postprocessing iter=" << iter << " violations=" << cur_violations;
    if (cur_violations == 0) break;
    if (cur_violations >= prev_violations) break; // no improvement, stop
    prev_violations = cur_violations;
  }

  // No fixing of nodes: postprocessing performs local moves to resolve
  // negative constraints and we no longer invoke the rebalancer/refinement
  // path here. Avoid manipulating FixedVertexSupport to prevent lifecycle
  // issues and potential null-pointer crashes.
  const auto constraints = readConstraintFile(context.partition.constraint_filename);

  LOG << "-------------- stats after postprocessing --------------";
  LOG << "km1       =" << metrics::quality(partitioned_hg, context);
  LOG << "Imbalance =" << metrics::imbalance(partitioned_hg, context);
  LOG << (verifyConstraints(partitioned_hg, context)
    ? "Constraints respected after postprocessing"
    : "!!! Constraints violated after postprocessing !!!");
  LOG << "";

  GainCachePtr::deleteGainCache(gain_cache);
}

} // namespace mt_kahypar::constraints