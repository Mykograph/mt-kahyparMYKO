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

#pragma once

#include <fstream>

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/factories.h"
#include "mt-kahypar/partition/refinement/gains/gain_cache_ptr.h"
#include "mt-kahypar/partition/refinement/i_rebalancer.h"
#include "mt-kahypar/datastructures/priority_queue.h"
#include "mt-kahypar/utils/cast.h"
#include "mt-kahypar/utils/exception.h"

namespace mt_kahypar::constraints {

using Key = std::tuple<HypernodeID, HypernodeID, HypernodeID>;

template<typename Comparator = std::less<Key>, uint32_t arity = 4>
using PQ = ds::Heap<Key, HypernodeID, Comparator, arity>;

// ----------------------------------------------------------------------------
// ConstraintGraph – flat adjacency list of anti‑affinity pairs
// ----------------------------------------------------------------------------
struct ConstraintGraph {
  explicit ConstraintGraph(const HypernodeID num_hypernodes) :
    adjacency(num_hypernodes) { }

  void addConstraint(const HypernodeID u, const HypernodeID v) {
    adjacency[u].push_back(v);
    adjacency[v].push_back(u);
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  const vec<HypernodeID>& neighbors(const HypernodeID u) const {
    return adjacency[u];
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HypernodeID degree(const HypernodeID u) const {
    return static_cast<HypernodeID>(adjacency[u].size());
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  bool isConstrained(const HypernodeID u) const {
    return !adjacency[u].empty();
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HypernodeID numNodes() const {
    return static_cast<HypernodeID>(adjacency.size());
  }

  vec<vec<HypernodeID>> adjacency;
};

namespace {

// Same line format as the constraint-file reader in hypergraph_factory.cc:
// whitespace separated "u v" pairs, one constraint per line.
inline void readConstraintPairs(const std::string& filename,
                                vec<std::pair<HypernodeID, HypernodeID>>& constraints) {
  std::ifstream file(filename);
  if ( !file.is_open() ) {
    throw InvalidInputException("Could not open constraint file: " + filename);
  }
  HypernodeID u = 0;
  HypernodeID v = 0;
  while ( file >> u >> v ) {
    constraints.emplace_back(u, v);
  }
  if ( file.bad() ) {
    throw InvalidInputException("Error while reading constraint file: " + filename);
  }
}

} // namespace

inline ConstraintGraph buildConstraintGraph(const std::string& constraint_file_name,
                                            const HypernodeID num_hypernodes) {
  vec<std::pair<HypernodeID, HypernodeID>> constraint_pairs;
  readConstraintPairs(constraint_file_name, constraint_pairs);

  ConstraintGraph graph(num_hypernodes);
  for ( const auto& [u, v] : constraint_pairs ) {
    if ( u >= num_hypernodes || v >= num_hypernodes ) {
      throw InvalidInputException(
        "Constraint file references node id out of range: (" + STR(u) + ", " + STR(v) + ")");
    }
    graph.addConstraint(u, v);
  }
  return graph;
}

// ----------------------------------------------------------------------------
// Helper: count violated constraints (pairs in the same block)
// ----------------------------------------------------------------------------
template<typename PartitionedHypergraph>
HypernodeID countViolatedConstraints(const PartitionedHypergraph& partitioned_hg,
                                     const ConstraintGraph& constraint_graph) {
  HypernodeID violated = 0;
  for (HypernodeID u = 0; u < constraint_graph.numNodes(); ++u) {
    if (!constraint_graph.isConstrained(u)) continue;
    const PartitionID part_u = partitioned_hg.partID(u);
    for (HypernodeID v : constraint_graph.neighbors(u)) {
      if (v <= u) continue;   // each undirected pair appears twice – count once
      if (partitioned_hg.partID(v) == part_u) {
        ++violated;
      }
    }
  }
  return violated;
}

template<typename PartitionedHypergraph>
bool verifyConstraints(const PartitionedHypergraph& partitioned_hg, const ConstraintGraph& constraint_graph) {
  for ( HypernodeID u = 0; u < constraint_graph.numNodes(); ++u ) {
    if ( !constraint_graph.isConstrained(u) ) continue;
    for ( const HypernodeID v : constraint_graph.neighbors(u) ) {
      if ( v <= u ) continue; // each undirected constraint is stored twice, check once
      if ( partitioned_hg.partID(u) == partitioned_hg.partID(v) ) {
        return false;
      }
    }
  }
  return true;
}

template<typename PartitionedHypergraph>
HypernodeID incidentNodesInSamePart(const PartitionedHypergraph& partitioned_hg,
                                    const ConstraintGraph& constraint_graph,
                                    const HypernodeID node_id) {
  HypernodeID num_nodes = 0;
  const PartitionID part = partitioned_hg.partID(node_id);
  for ( const HypernodeID neighbor : constraint_graph.neighbors(node_id) ) {
    if ( partitioned_hg.partID(neighbor) == part ) {
      ++num_nodes;
    }
  }
  return num_nodes;
}

template<typename PartitionedHypergraph>
PartitionID getBestCutPartition(const bool must_cut_be_positive,
                                const HypernodeID& node_id,
                                const PartitionID& current_partition,
                                const vec<bool>& is_partition_invalid,
                                gain_cache_t& gain_cache) {
  PartitionID best_partition = current_partition;
  const PartitionID num_partitons = is_partition_invalid.size();
  HyperedgeWeight max_gain = std::numeric_limits<HyperedgeWeight>::min();

  for (PartitionID partition = 0; partition < num_partitons; partition++) {
    if (is_partition_invalid[partition]) continue;
    HyperedgeWeight gain = GainCachePtr::applyWithConcreteGainCacheForHG<PartitionedHypergraph>(
      [&](const auto& cache){
        return cache.gain(node_id, current_partition, partition);
      },
      gain_cache
    );
    if (must_cut_be_positive && gain < 0) continue;
    if (gain > max_gain) {
      max_gain = gain;
      best_partition = partition;
    }
  }
  return best_partition;
}

template<typename PartitionedHypergraph>
void descendingConstraintDegree(PartitionedHypergraph& partitioned_hg,
                                const Context& context,
                                const ConstraintGraph& constraint_graph,
                                gain_cache_t& gain_cache) {
  unused(context);
  GainCachePtr::applyWithConcreteGainCacheForHG<PartitionedHypergraph>(
    [&](auto& cache){
        cache.initializeGainCache(partitioned_hg);
    },
    gain_cache
  );
  auto delta_func =
  [&partitioned_hg, &gain_cache](const SynchronizedEdgeUpdate& sync_update) {
    GainCachePtr::applyWithConcreteGainCacheForHG<PartitionedHypergraph>(
      [&](auto& cache) {
        cache.deltaGainUpdate(partitioned_hg, sync_update);
      },
      gain_cache
    );
  };

  // initialize PQ -- only constrained nodes are inserted
  std::vector<PosT> positions(constraint_graph.numNodes(), invalid_position);
  PQ<std::less<Key>, 4> heap(positions.data(), positions.size());
  for ( HypernodeID node = 0; node < constraint_graph.numNodes(); ++node ) {
    if ( !constraint_graph.isConstrained(node) ) continue;
    heap.insert(node, { incidentNodesInSamePart(partitioned_hg, constraint_graph, node),
                        constraint_graph.degree(node), node });
  }

  while ( !heap.empty() ) {
    const HypernodeID node_id = heap.top();
    heap.deleteTop();
    const PartitionID partition_id = partitioned_hg.partID(node_id);
    vec<bool> invalid_partitions(partitioned_hg.k(), false);

    for ( const HypernodeID neighbor : constraint_graph.neighbors(node_id) ) {
      const PartitionID neighbor_partition_id = partitioned_hg.partID(neighbor);
      invalid_partitions[neighbor_partition_id] = true;
    }
    const PartitionID new_partition_id = getBestCutPartition<PartitionedHypergraph>(
      false, node_id, partition_id, invalid_partitions, gain_cache);
    if ( new_partition_id != partition_id ) {
      partitioned_hg.changeNodePart(node_id, partition_id, new_partition_id, delta_func);
    }
  }
}

template<typename PartitionedHypergraph>
void postprocessNegativeConstraints(PartitionedHypergraph& partitioned_hg,
                                    const Context& context) {
  if ( context.partition.constraint_file_name.empty() ) {
    return;
  }

  ConstraintGraph constraint_graph = buildConstraintGraph(
    context.partition.constraint_file_name, partitioned_hg.initialNumNodes());

  gain_cache_t gain_cache = GainCachePtr::constructGainCache(context);
  std::unique_ptr<IRebalancer> rebalancer = RebalancerFactory::getInstance().createObject(
      context.refinement.rebalancing.algorithm, partitioned_hg.initialNumNodes(), context, gain_cache);

  // ---- initial violation count ----
  HypernodeID initial_violations = countViolatedConstraints(partitioned_hg, constraint_graph);
  LOG << "Initial violated constraints before any fix: " << initial_violations;

  // ---- iterative repair loop ----
  const size_t max_iterations = 50;   // safety guard
  size_t iteration = 0;
  bool constraints_ok = false;
  bool balance_ok = false;

  while (iteration < max_iterations) {
    ++iteration;
    LOG << "Iteration " << iteration << ": fixing constraints...";

    // 1. Fix constraints (descending constraint degree)
    descendingConstraintDegree(partitioned_hg, context, constraint_graph, gain_cache);

    // 2. Check constraints and balance after fixing
    HypernodeID violations = countViolatedConstraints(partitioned_hg, constraint_graph);
    double imbalance = metrics::imbalance(partitioned_hg, context);
    LOG << "After fix: violations = " << violations << ", imbalance = " << imbalance;

    constraints_ok = (violations == 0);
    balance_ok = (imbalance <= context.partition.epsilon + 1e-9);

    if (constraints_ok && balance_ok) {
      LOG << "Both constraints and balance satisfied after fix.";
      break;
    }

    // 3. If constraints are ok but balance is not, run rebalancer
    if (constraints_ok && !balance_ok) {
      LOG << "Constraints ok, running rebalancer to fix imbalance...";
      Metrics metrics { metrics::quality(partitioned_hg, context), imbalance };
      mt_kahypar_partitioned_hypergraph_t phg = utils::partitioned_hg_cast(partitioned_hg);
      rebalancer->initialize(phg);
      rebalancer->refine(phg, {}, metrics, 0.0);

      // 4. Re-check after rebalancer
      violations = countViolatedConstraints(partitioned_hg, constraint_graph);
      imbalance = metrics::imbalance(partitioned_hg, context);
      LOG << "After rebalancer: violations = " << violations << ", imbalance = " << imbalance;

      constraints_ok = (violations == 0);
      balance_ok = (imbalance <= context.partition.epsilon + 1e-9);

      if (constraints_ok && balance_ok) {
        LOG << "Both constraints and balance satisfied after rebalancer.";
        break;
      }
    } else {
      // If constraints are not ok, rebalancing might break them further, so we loop again.
      // But we can still try rebalancing if imbalance is bad, but better to loop to fix constraints first.
      // However, sometimes rebalancing helps with constraints indirectly, so we can run it anyway.
      LOG << "Constraints not ok, rebalancing might help or hurt; will try rebalancer anyway.";
      Metrics metrics { metrics::quality(partitioned_hg, context), imbalance };
      mt_kahypar_partitioned_hypergraph_t phg = utils::partitioned_hg_cast(partitioned_hg);
      rebalancer->initialize(phg);
      rebalancer->refine(phg, {}, metrics, 0.0);

      // After rebalancer, re-check (it may have improved or worsened constraints)
      violations = countViolatedConstraints(partitioned_hg, constraint_graph);
      imbalance = metrics::imbalance(partitioned_hg, context);
      LOG << "After rebalancer: violations = " << violations << ", imbalance = " << imbalance;

      constraints_ok = (violations == 0);
      balance_ok = (imbalance <= context.partition.epsilon + 1e-9);

      if (constraints_ok && balance_ok) {
        LOG << "Both constraints and balance satisfied after rebalancer.";
        break;
      }
    }
  }

  // ---- final log ----
  LOG << "-------------- final stats after postprocessing --------------";
  LOG << "km1       = " << metrics::quality(partitioned_hg, context);
  LOG << "Imbalance = " << metrics::imbalance(partitioned_hg, context);
  LOG << "Constraints satisfied: " << (verifyConstraints(partitioned_hg, constraint_graph) ? "yes" : "no");
  if (!constraints_ok || !balance_ok) {
    LOG << "WARNING: Not all constraints satisfied or imbalance exceeds epsilon after " << iteration << " iterations!";
  } else {
    LOG << "SUCCESS: All constraints satisfied and imbalance within epsilon.";
  }
  LOG << "";

  GainCachePtr::deleteGainCache(gain_cache);
}

} // namespace mt_kahypar::constraints