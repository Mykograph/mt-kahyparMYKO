#pragma once

#include <fstream>
#include <limits>
#include <memory>
#include <vector>
#include <algorithm>
#include <deque>

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/factories.h"
#include "mt-kahypar/partition/refinement/gains/gain_cache_ptr.h"
#include "mt-kahypar/datastructures/priority_queue.h"
#include "mt-kahypar/utils/cast.h"

namespace mt_kahypar::constraints {

// ####################### Helpers #######################

// Reads constraint pairs directly from file — no dependency on io:: functions.
//
// max_valid_id / num_skipped let callers detect out-of-range IDs, e.g. a
// generator that emitted 1-indexed .hgr vertex numbers instead of
// mt-kahypar's 0-indexed internal HypernodeIDs. Previously any u/v >= 0 was
// accepted unconditionally and used directly to index into partID() — an
// out-of-range ID silently read garbage instead of failing loudly.
inline vec<std::pair<HypernodeID, HypernodeID>> readConstraintFile(
    const std::string& filename,
    HypernodeID max_valid_id = std::numeric_limits<HypernodeID>::max(),
    size_t* num_skipped = nullptr) {
  vec<std::pair<HypernodeID, HypernodeID>> constraints;
  if (num_skipped) *num_skipped = 0;
  if (filename.empty()) return constraints;
  std::ifstream file(filename);
  if (!file) return constraints;
  int64_t u = 0, v = 0;
  size_t skipped = 0;
  while (file >> u >> v) {
    if (u >= 0 && v >= 0 && u != v) {
      if (static_cast<uint64_t>(u) > static_cast<uint64_t>(max_valid_id) ||
          static_cast<uint64_t>(v) > static_cast<uint64_t>(max_valid_id)) {
        ++skipped;
        continue;
      }
      constraints.emplace_back(
        static_cast<HypernodeID>(u),
        static_cast<HypernodeID>(v));
    }
  }
  if (num_skipped) *num_skipped = skipped;
  return constraints;
}

// Convenience: derive max_valid_id from a partitioned hypergraph.
template<typename PartitionedHypergraph>
inline HypernodeID maxValidNodeId(const PartitionedHypergraph& partitioned_hg) {
  const size_t n = partitioned_hg.initialNumNodes();
  return n == 0 ? std::numeric_limits<HypernodeID>::max()
                : static_cast<HypernodeID>(n - 1);
}

// ####################### Verify #######################

template<typename PartitionedHypergraph>
bool verifyConstraints(const PartitionedHypergraph& partitioned_hg, const Context& context) {
  const auto constraints = readConstraintFile(
      context.partition.constraint_filename, maxValidNodeId(partitioned_hg));
  for (const auto& [u, v] : constraints) {
    if (partitioned_hg.partID(u) == partitioned_hg.partID(v)) {
      return false;
    }
  }
  return true;
}

template<typename PartitionedHypergraph>
size_t countConstraintViolations(const PartitionedHypergraph& partitioned_hg, const Context& context) {
  const auto constraints = readConstraintFile(
      context.partition.constraint_filename, maxValidNodeId(partitioned_hg));
  size_t cnt = 0;
  for (const auto& [u, v] : constraints) {
    if (partitioned_hg.partID(u) == partitioned_hg.partID(v)) ++cnt;
  }
  return cnt;
}

// ####################### Diagnostics #######################

// Connected-component size histogram for the constraint graph. A component
// of size <= k is always resolvable to zero violations regardless of its
// internal structure. A component of size > k is *not* automatically
// infeasible (e.g. a large bipartite component is still 2-colorable) but is
// flagged here as worth inspecting, since a fully-connected component
// (clique) of that size would need exactly that many colors to be
// violation-free at all.
inline void logConstraintGraphComponentSizes(const ds::DynamicGraph& cg, PartitionID k) {
  const size_t n = cg.initialNumNodes();
  vec<char> visited(n, 0);
  size_t num_components = 0;
  size_t max_component = 0;
  size_t components_over_k = 0;
  std::deque<HypernodeID> bfs;

  for (HypernodeID start : cg.nodes()) {
    if (visited[start]) continue;
    ++num_components;
    size_t size = 0;
    bfs.push_back(start);
    visited[start] = 1;
    while (!bfs.empty()) {
      const HypernodeID cur = bfs.front();
      bfs.pop_front();
      ++size;
      for (const auto& e : cg.incidentEdges(cur)) {
        const HypernodeID nb = cg.edge(e).target;
        if (!visited[nb]) {
          visited[nb] = 1;
          bfs.push_back(nb);
        }
      }
    }
    max_component = std::max(max_component, size);
    if (size > static_cast<size_t>(k)) ++components_over_k;
  }

  LOG << "Constraint graph components=" << num_components
      << " max_component_size=" << max_component
      << " components_larger_than_k(k=" << k << ")=" << components_over_k;
  if (components_over_k > 0) {
    LOG << "NOTE: components larger than k are not guaranteed violation-free. "
        << "If such a component is (close to) a clique of size m, the minimum "
        << "possible violated pairs at this k is C(m,2) minus the sum of "
        << "C(group_size,2) over the most even k-way split of that clique — "
        << "i.e. some violations there may be mathematically unavoidable, not a bug.";
  }
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

// Propagation-based fixer: process a queue of violating constraint nodes and
// move them to valid partitions (not occupied by any constraint neighbor).
// This is incremental and reacts to neighbor moves until no violations
// remain or a move limit is reached.
template<typename PartitionedHypergraph>
void propagateConstraintFixes(PartitionedHypergraph& partitioned_hg,
                             gain_cache_t& gain_cache,
                             const Context& context) {
  GainCachePtr::applyWithConcreteGainCacheForHG<PartitionedHypergraph>(
    [&](auto& cache) { cache.initializeGainCache(partitioned_hg); }, gain_cache);

  auto delta_func = [&partitioned_hg, &gain_cache](const SynchronizedEdgeUpdate& sync_update) {
    GainCachePtr::applyWithConcreteGainCacheForHG<PartitionedHypergraph>([&](auto& cache) {
      cache.deltaGainUpdate(partitioned_hg, sync_update);
    }, gain_cache);
  };

  const ds::DynamicGraph& cg = const_cast<const PartitionedHypergraph&>(partitioned_hg).fixedVertexSupport().getConstraintGraph();
  const size_t n = cg.initialNumNodes();
  std::deque<HypernodeID> q;
  vec<char> in_queue(n, 0);

  // Initialize queue with all currently violating constraint nodes
  for (HypernodeID cnode : cg.nodes()) {
    if (incidentNodesInSamePart(partitioned_hg, cnode) > 0) {
      q.push_back(cnode);
      in_queue[cnode] = 1;
    }
  }

  const size_t max_moves = std::max<size_t>(n * 4, 1000);
  size_t moves = 0;

  while (!q.empty() && moves < max_moves) {
    HypernodeID cnode = q.front(); q.pop_front(); in_queue[cnode] = 0;
    const HypernodeID hg_node_id = static_cast<HypernodeID>(cg.nodeWeight(cnode));
    const PartitionID cur_part = partitioned_hg.partID(hg_node_id);

    // compute invalid partitions occupied by neighbors
    vec<bool> invalid_partitions(partitioned_hg.k(), false);
    for ( const auto& e : cg.incidentEdges(cnode) ) {
      HypernodeID nb = cg.edge(e).target;
      HypernodeID nb_hg = static_cast<HypernodeID>(cg.nodeWeight(nb));
      invalid_partitions[ partitioned_hg.partID(nb_hg) ] = true;
    }

    // still violating?
    if (!invalid_partitions[cur_part]) continue;

    // choose best partition avoiding neighbors' blocks
    PartitionID best = getBestPartition<PartitionedHypergraph>(hg_node_id, cur_part, invalid_partitions, gain_cache);
    if (best == cur_part) continue;

    partitioned_hg.changeNodePart(hg_node_id, cur_part, best, delta_func);
    ++moves;

    // after a move, enqueue all constraint-graph neighbors (they may now violate)
    for ( const auto& e : cg.incidentEdges(cnode) ) {
      HypernodeID nb = cg.edge(e).target;
      if (!in_queue[nb] && incidentNodesInSamePart(partitioned_hg, nb) > 0) {
        q.push_back(nb);
        in_queue[nb] = 1;
      }
    }

    // also re-enqueue the moved node if it still violates
    if (incidentNodesInSamePart(partitioned_hg, cnode) > 0 && !in_queue[cnode]) {
      q.push_back(cnode);
      in_queue[cnode] = 1;
    }
  }

  LOG << "propagateConstraintFixes moved=" << moves << " max_moves=" << max_moves;
}

// ####################### Entry point #######################

template<typename PartitionedHypergraph>
void postprocessNegativeConstraints(PartitionedHypergraph& partitioned_hg,
                                    const Context& context) {
  gain_cache_t gain_cache = GainCachePtr::constructGainCache(context);
  LOG << "Starting postprocessNegativeConstraints; constraint_file='" << context.partition.constraint_filename << "'";

  // ---- Diagnostics: run before any fixing, so results reflect the raw input ----
  {
    size_t skipped_oob = 0;
    const auto constraints_check = readConstraintFile(
        context.partition.constraint_filename, maxValidNodeId(partitioned_hg), &skipped_oob);
    LOG << "Constraint file entries kept=" << constraints_check.size()
        << " | skipped out-of-range=" << skipped_oob
        << " | valid HypernodeID range=[0, " << maxValidNodeId(partitioned_hg) << "]";
    if (skipped_oob > 0) {
      LOG << "!!! " << skipped_oob << " constraint entries referenced HypernodeIDs outside the "
          << "valid range and were skipped. This usually means the constraint generator used "
          << "1-indexed IDs (matching the .hgr file's own numbering) instead of mt-kahypar's "
          << "0-indexed internal HypernodeIDs. Check the generator before trusting any results below.";
    }
  }

  try {
    const ds::DynamicGraph& dbg_cg = const_cast<const PartitionedHypergraph&>(partitioned_hg).fixedVertexSupport().getConstraintGraph();
    size_t num_cnodes = dbg_cg.initialNumNodes();
    size_t init_viol = 0;
    for (HypernodeID cn : dbg_cg.nodes()) if (incidentNodesInSamePart(partitioned_hg, cn) > 0) ++init_viol;
    LOG << "Constraint graph nodes=" << num_cnodes << " initial_violations=" << init_viol;
    logConstraintGraphComponentSizes(dbg_cg, partitioned_hg.k());
  } catch (...) {
    LOG << "No constraint graph attached or error while inspecting it.";
  }

  // Use a propagation-based fixer that incrementally processes violating
  // constraint nodes and reacts to neighbor moves until convergence.
  // If the simple propagation approach fails, use a more robust iterative
  // fixer that pins moved nodes temporarily and enforces block capacity
  // constraints to avoid partition collapse.
  auto robustConstraintFixes = [&] (PartitionedHypergraph& phg, gain_cache_t& gc, const Context& ctx) {
    GainCachePtr::applyWithConcreteGainCacheForHG<PartitionedHypergraph>(
      [&](auto& cache) { cache.initializeGainCache(phg); }, gc);

    auto delta_func_local = [&phg, &gc](const SynchronizedEdgeUpdate& sync_update) {
      GainCachePtr::applyWithConcreteGainCacheForHG<PartitionedHypergraph>([&](auto& cache) {
        cache.deltaGainUpdate(phg, sync_update);
      }, gc);
    };

    const ds::DynamicGraph& cg = const_cast<const PartitionedHypergraph&>(phg).fixedVertexSupport().getConstraintGraph();
    const size_t n = cg.initialNumNodes();

    // Raised from 16: the previous hard cap could silently stop the fixer
    // while it was still making progress. We now also log whether we exit
    // via convergence or via hitting this cap, so that distinction is no
    // longer invisible in the output.
    const size_t max_rounds = 200;
    const size_t pin_rounds = 2; // number of rounds to pin moved nodes
    vec<uint8_t> pin_cooldown(n, 0);

    size_t total_moves = 0;
    size_t round = 0;
    for (; round < max_rounds; ++round) {
      // collect violating nodes
      vec<std::pair<HypernodeID, HypernodeID>> viols; // (cnode, violation_count)
      for (HypernodeID cnode : cg.nodes()) {
        if (incidentNodesInSamePart(phg, cnode) > 0) {
          viols.emplace_back(cnode, incidentNodesInSamePart(phg, cnode));
        }
      }

      LOG << "  [robustConstraintFixes] round=" << round << " violating_nodes=" << viols.size();
      if (viols.empty()) break;

      // process highest violation nodes first
      std::sort(viols.begin(), viols.end(), [](const auto& a, const auto& b){ return a.second > b.second; });

      size_t moves_this_round = 0;
      for (const auto& [cnode, cnt] : viols) {
        if (pin_cooldown[cnode] > 0) continue;
        const HypernodeID hg_node_id = static_cast<HypernodeID>(cg.nodeWeight(cnode));
        const PartitionID cur_part = phg.partID(hg_node_id);

        // compute invalid partitions
        vec<bool> invalid_partitions(phg.k(), false);
        for (const auto& e : cg.incidentEdges(cnode)) {
          HypernodeID nb = cg.edge(e).target;
          HypernodeID nb_hg = static_cast<HypernodeID>(cg.nodeWeight(nb));
          invalid_partitions[ phg.partID(nb_hg) ] = true;
        }
        if (!invalid_partitions[cur_part]) continue; // no longer violating

        // evaluate best allowed partition under capacity
        PartitionID best = cur_part;
        HyperedgeWeight best_gain = std::numeric_limits<HyperedgeWeight>::min();
        for (PartitionID p = 0; p < phg.k(); ++p) {
          if (invalid_partitions[p]) continue;
          const HypernodeWeight node_w = phg.nodeWeight(hg_node_id);
          if (phg.partWeight(p) + node_w > ctx.partition.max_part_weights[p]) continue;
          HyperedgeWeight gain = GainCachePtr::applyWithConcreteGainCacheForHG<PartitionedHypergraph>(
            [&](const auto& cache) {
              return cache.gain(hg_node_id, cur_part, p);
            }, gc);
          if (gain > best_gain) { best_gain = gain; best = p; }
        }
        if (best == cur_part) continue;

        // attempt move with gain-cache-aware changeNodePart
        bool moved = phg.changeNodePart(hg_node_id, cur_part, best,
          ctx.partition.max_part_weights[best], [&]() { ++moves_this_round; ++total_moves; }, delta_func_local);
        if (moved) {
          pin_cooldown[cnode] = pin_rounds;
          // also lightly pin neighbors to avoid immediate oscillation
          for (const auto& e : cg.incidentEdges(cnode)) {
            HypernodeID nb = cg.edge(e).target;
            pin_cooldown[nb] = std::max<uint8_t>(pin_cooldown[nb], 1);
          }
        }
      }

      // decrease cooldowns
      for (size_t i = 0; i < n; ++i) if (pin_cooldown[i] > 0) --pin_cooldown[i];

      if (moves_this_round == 0) break;
    }

    if (round == max_rounds) {
      LOG << "robustConstraintFixes: hit max_rounds=" << max_rounds
          << " while violations were still present — this run may not have converged; "
          << "consider raising max_rounds further.";
    }
    LOG << "robustConstraintFixes total_moves=" << total_moves << " rounds_used=" << (round + 1);
  };

  robustConstraintFixes(partitioned_hg, gain_cache, context);

  // No fixing of nodes: postprocessing performs local moves to resolve
  // negative constraints and we no longer invoke the rebalancer/refinement
  // path here. Avoid manipulating FixedVertexSupport to prevent lifecycle
  // issues and potential null-pointer crashes.

  LOG << "-------------- stats after postprocessing --------------";
  LOG << "km1       =" << metrics::quality(partitioned_hg, context);
  LOG << "Imbalance =" << metrics::imbalance(partitioned_hg, context);
  LOG << (verifyConstraints(partitioned_hg, context)
    ? "Constraints respected after postprocessing"
    : "!!! Constraints violated after postprocessing !!!");
  LOG << "Remaining violated constraints=" << countConstraintViolations(partitioned_hg, context);
  LOG << "";

  GainCachePtr::deleteGainCache(gain_cache);
}

} // namespace mt_kahypar::constraints