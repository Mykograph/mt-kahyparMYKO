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
#include <deque>
#include <algorithm>
#include <limits>
#include <stdexcept>

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/factories.h"
#include "mt-kahypar/partition/refinement/gains/gain_cache_ptr.h"
#include "mt-kahypar/partition/refinement/i_rebalancer.h"
#include "mt-kahypar/datastructures/priority_queue.h"
#include "mt-kahypar/utils/cast.h"
#include "mt-kahypar/utils/exception.h"

namespace mt_kahypar::constraints {

// ConstraintGraph – flat adjacency list 
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

// Reads constraint pairs from file 
inline void readConstraintPairs(const std::string& filename,
                                vec<std::pair<HypernodeID, HypernodeID>>& constraints) {
  std::ifstream file(filename);
  if ( !file.is_open() ) {
    throw InvalidInputException("Could not open constraint file: " + filename);
  }
  HypernodeID u = 0, v = 0;
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

// Helper: count violated constraints 
template<typename PartitionedHypergraph>
HypernodeID countViolatedConstraints(const PartitionedHypergraph& partitioned_hg,
                                     const ConstraintGraph& constraint_graph) {
  HypernodeID violated = 0;
  for (HypernodeID u = 0; u < constraint_graph.numNodes(); ++u) {
    if (!constraint_graph.isConstrained(u)) continue;
    const PartitionID part_u = partitioned_hg.partID(u);
    for (HypernodeID v : constraint_graph.neighbors(u)) {
      if (v <= u) continue;   // count each undirected pair once
      if (partitioned_hg.partID(v) == part_u) {
        ++violated;
      }
    }
  }
  return violated;
}

template<typename PartitionedHypergraph>
bool verifyConstraints(const PartitionedHypergraph& partitioned_hg,
                       const ConstraintGraph& constraint_graph) {
  for (HypernodeID u = 0; u < constraint_graph.numNodes(); ++u) {
    if (!constraint_graph.isConstrained(u)) continue;
    for (const HypernodeID v : constraint_graph.neighbors(u)) {
      if (v <= u) continue;
      if (partitioned_hg.partID(u) == partitioned_hg.partID(v)) {
        return false;
      }
    }
  }
  return true;
}

// Core fixing helpers 

// Counts how many constraint neighbours of a node are in the same block
template<typename PartitionedHypergraph>
HypernodeID incidentNodesInSamePart(const PartitionedHypergraph& partitioned_hg,
                                    const ConstraintGraph& constraint_graph,
                                    const HypernodeID node_id) {
  HypernodeID num = 0;
  const PartitionID part = partitioned_hg.partID(node_id);
  for (const HypernodeID neighbor : constraint_graph.neighbors(node_id)) {
    if (partitioned_hg.partID(neighbor) == part) ++num;
  }
  return num;
}

// Picks the partition with maximum gain, respecting:
//   - invalid_partitions 
//   - capacity limits 
template<typename PartitionedHypergraph>
PartitionID getBestPartitionWithCapacity(const HypernodeID node_id,
                                         const PartitionID current_partition,
                                         const vec<bool>& invalid_partitions,
                                         const PartitionedHypergraph& partitioned_hg,
                                         const Context& context,
                                         gain_cache_t& gain_cache) {
  PartitionID best = current_partition;
  HyperedgeWeight best_gain = std::numeric_limits<HyperedgeWeight>::min();
  const HypernodeWeight node_w = partitioned_hg.nodeWeight(node_id);

  for (PartitionID p = 0; p < partitioned_hg.k(); ++p) {
    if (invalid_partitions[p]) continue;
    // Capacity check – only allow move if target partition can accommodate the node
    if (partitioned_hg.partWeight(p) + node_w > context.partition.max_part_weights[p])
      continue;

    const HyperedgeWeight gain = GainCachePtr::applyWithConcreteGainCacheForHG<PartitionedHypergraph>(
      [&](const auto& cache) { return cache.gain(node_id, current_partition, p); },
      gain_cache
    );
    if (gain > best_gain) {
      best_gain = gain;
      best = p;
    }
  }
  return best;
}


// Propagation-based fixer – processes a queue of violating nodes and
// propagates moves to neighbours.
template<typename PartitionedHypergraph>
void propagateConstraintFixes(PartitionedHypergraph& partitioned_hg,
                              ConstraintGraph& constraint_graph,
                              gain_cache_t& gain_cache,
                              const Context& context) {
  // Initialise gain cache
  GainCachePtr::applyWithConcreteGainCacheForHG<PartitionedHypergraph>(
    [&](auto& cache) { cache.initializeGainCache(partitioned_hg); },
    gain_cache
  );

  auto delta_func = [&partitioned_hg, &gain_cache](const SynchronizedEdgeUpdate& sync_update) {
    GainCachePtr::applyWithConcreteGainCacheForHG<PartitionedHypergraph>(
      [&](auto& cache) { cache.deltaGainUpdate(partitioned_hg, sync_update); },
      gain_cache
    );
  };

  const size_t n = constraint_graph.numNodes();
  std::deque<HypernodeID> q;
  vec<char> in_queue(n, 0);

  // Initialise queue with all violating nodes
  for (HypernodeID node = 0; node < n; ++node) {
    if (constraint_graph.isConstrained(node) &&
        incidentNodesInSamePart(partitioned_hg, constraint_graph, node) > 0) {
      q.push_back(node);
      in_queue[node] = 1;
    }
  }

  const size_t max_moves = std::max<size_t>(n * 4, 1000);
  size_t moves = 0;

  while (!q.empty() && moves < max_moves) {
    HypernodeID node = q.front(); q.pop_front(); in_queue[node] = 0;

    const PartitionID cur_part = partitioned_hg.partID(node);
    // Determine partitions occupied by constraint neighbours
    vec<bool> invalid(partitioned_hg.k(), false);
    for (const HypernodeID nb : constraint_graph.neighbors(node)) {
      invalid[partitioned_hg.partID(nb)] = true;
    }

    // If this node is no longer violating, skip
    if (!invalid[cur_part]) continue;

    // Choose best valid partition (respecting capacities)
    const PartitionID best = getBestPartitionWithCapacity(
        node, cur_part, invalid, partitioned_hg, context, gain_cache);

    if (best == cur_part) continue;  // no move possible

    // Perform the move
    partitioned_hg.changeNodePart(node, cur_part, best, delta_func);
    ++moves;

    // Enqueue all neighbours (they might now be in conflict)
    for (const HypernodeID nb : constraint_graph.neighbors(node)) {
      if (!in_queue[nb] &&
          incidentNodesInSamePart(partitioned_hg, constraint_graph, nb) > 0) {
        q.push_back(nb);
        in_queue[nb] = 1;
      }
    }

    // Re‑enqueue the moved node if it still violates
    if (incidentNodesInSamePart(partitioned_hg, constraint_graph, node) > 0 && !in_queue[node]) {
      q.push_back(node);
      in_queue[node] = 1;
    }
  }

  LOG << "propagateConstraintFixes moved=" << moves << " max_moves=" << max_moves;
}

// Robust iterative fixer with pinning 
template<typename PartitionedHypergraph>
void robustConstraintFixes(PartitionedHypergraph& partitioned_hg,
                           ConstraintGraph& constraint_graph,
                           gain_cache_t& gain_cache,
                           const Context& context) {
  GainCachePtr::applyWithConcreteGainCacheForHG<PartitionedHypergraph>(
    [&](auto& cache) { cache.initializeGainCache(partitioned_hg); },
    gain_cache
  );

  auto delta_func = [&partitioned_hg, &gain_cache](const SynchronizedEdgeUpdate& sync_update) {
    GainCachePtr::applyWithConcreteGainCacheForHG<PartitionedHypergraph>(
      [&](auto& cache) { cache.deltaGainUpdate(partitioned_hg, sync_update); },
      gain_cache
    );
  };

  const size_t n = constraint_graph.numNodes();
  const size_t max_rounds = 200;
  const size_t pin_rounds = 2;   // how many rounds a moved node is frozen
  vec<uint8_t> pin_cooldown(n, 0);

  size_t total_moves = 0;
  size_t round = 0;

  for (; round < max_rounds; ++round) {
    // Collect violating nodes with their violation counts
    vec<std::pair<HypernodeID, HypernodeID>> viols;
    for (HypernodeID node = 0; node < n; ++node) {
      if (constraint_graph.isConstrained(node)) {
        const HypernodeID cnt = incidentNodesInSamePart(partitioned_hg, constraint_graph, node);
        if (cnt > 0) viols.emplace_back(node, cnt);
      }
    }
    if (viols.empty()) break;

    // Process worst offenders first
    std::sort(viols.begin(), viols.end(),
              [](const auto& a, const auto& b) { return a.second > b.second; });

    size_t moves_this_round = 0;
    for (const auto& [node, cnt] : viols) {
      if (pin_cooldown[node] > 0) continue;

      const PartitionID cur_part = partitioned_hg.partID(node);

      vec<bool> invalid(partitioned_hg.k(), false);
      for (const HypernodeID nb : constraint_graph.neighbors(node)) {
        invalid[partitioned_hg.partID(nb)] = true;
      }

      if (!invalid[cur_part]) continue;   // not actually violating now (may have been fixed)

      const PartitionID best = getBestPartitionWithCapacity(
          node, cur_part, invalid, partitioned_hg, context, gain_cache);

      if (best == cur_part) continue;

      partitioned_hg.changeNodePart(node, cur_part, best, delta_func);
      ++moves_this_round;
      ++total_moves;

      // Freeze this node and its neighbours temporarily
      pin_cooldown[node] = pin_rounds;
      for (const HypernodeID nb : constraint_graph.neighbors(node)) {
        pin_cooldown[nb] = std::max<uint8_t>(pin_cooldown[nb], 1);
      }
    }

    // Decrease cooldowns
    for (size_t i = 0; i < n; ++i) {
      if (pin_cooldown[i] > 0) --pin_cooldown[i];
    }

    if (moves_this_round == 0) break;
  }

  LOG << "robustConstraintFixes total_moves=" << total_moves
      << " rounds_used=" << (round + 1);
  if (round == max_rounds) {
    LOG << "WARNING: reached max_rounds=" << max_rounds
        << " – consider raising further.";
  }
}

// Main postprocessing entry point 
template<typename PartitionedHypergraph>
void postprocessNegativeConstraints(PartitionedHypergraph& partitioned_hg,
                                    Context& context) {
  if (context.partition.constraint_file_name.empty()) {
    return;
  }

  // Build constraint graph from file
  ConstraintGraph constraint_graph =
      buildConstraintGraph(context.partition.constraint_file_name,
                           partitioned_hg.initialNumNodes());

  // Log initial violations
  const HypernodeID init_viol = countViolatedConstraints(partitioned_hg, constraint_graph);
  LOG << "Initial violated constraints: " << init_viol;
  context.partition.violated_constraints_after_refinement = init_viol;

  // Construct gain cache
  gain_cache_t gain_cache = GainCachePtr::constructGainCache(context);

  // First, try the fast propagation fixer
  propagateConstraintFixes(partitioned_hg, constraint_graph, gain_cache, context);

  // If violations remain, run the more robust iterative fixer
  const HypernodeID viol_after_prop = countViolatedConstraints(partitioned_hg, constraint_graph);
  if (viol_after_prop > 0) {
    LOG << "After propagation, still " << viol_after_prop
        << " violations – running robust fixer.";
    robustConstraintFixes(partitioned_hg, constraint_graph, gain_cache, context);
  }

  // Final verification
  const HypernodeID final_viol = countViolatedConstraints(partitioned_hg, constraint_graph);
  const double imbalance = metrics::imbalance(partitioned_hg, context);
  LOG << "-------------- final stats --------------";
  LOG << "km1       = " << metrics::quality(partitioned_hg, context);
  LOG << "Imbalance = " << imbalance;
  LOG << "Violations = " << final_viol;

  // Clean up gain cache 
  GainCachePtr::deleteGainCache(gain_cache);

  if (final_viol == 0 && imbalance <= context.partition.epsilon + 1e-9) {
    LOG << "SUCCESS: All constraints satisfied and balance within epsilon.";
  } else {
    std::string msg = "postprocessNegativeConstraints failed to satisfy all constraints: ";
    if (final_viol > 0) {
      msg += STR(final_viol) + " constraint violation(s) remain";
    }
    if (imbalance > context.partition.epsilon + 1e-9) {
      if (final_viol > 0) {
        msg += "; ";
      }
      msg += "imbalance " + STR(imbalance) + " exceeds epsilon " + STR(context.partition.epsilon);
    }
    LOG << "ERROR: " << msg;
    throw std::runtime_error(msg);
  }
}

} // namespace mt_kahypar::constraints