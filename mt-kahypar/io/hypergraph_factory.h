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
#include <string>
#include "include/mtkahypartypes.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/partition/context_enum_classes.h"
#include "mt-kahypar/utils/cast.h"
#include "mt-kahypar/io/hypergraph_io.h"  // for HyperedgeVector, EdgeVector typedefs used by OriginalEdgeSnapshot

namespace mt_kahypar {
namespace io {

// Snapshot of the edge list / weights exactly as parsed from the input file,
// taken BEFORE the constraint pass reweights existing edges or appends
// synthetic constraint edges. Pass a pointer to readInputFile(...) to have it
// populated, then use countOriginalCut(...) together with a final partition
// to compute the cut on the original, un-tampered instance.
struct OriginalEdgeSnapshot {
  bool is_graph = false;                   // true -> original_edges is populated, false -> original_hyperedges
  bool valid = false;                      // true once readInputFile has populated this snapshot
  HyperedgeVector original_hyperedges;     // populated when is_graph == false
  EdgeVector original_edges;                // populated when is_graph == true
  vec<HyperedgeWeight> original_edge_weights;
};

mt_kahypar_hypergraph_t readInputFile(const std::string& filename,
                                      const PresetType& preset,
                                      const InstanceType& instance,
                                      const FileFormat& format,
                                      const bool stable_construction,
                                      const bool remove_single_pin_hes,
                                      const bool print_warnings,
                                      const std::string& constraint_filename = "",
                                      const HyperedgeWeight constraint_weight = -100,
                                      const HyperedgeWeight normal_edge_weight_bonus = 0,
                                      OriginalEdgeSnapshot* out_snapshot = nullptr);

template<typename Hypergraph>
Hypergraph readInputFile(const std::string& filename,
                         const FileFormat& format,
                         const bool stable_construction,
                         const bool remove_single_pin_hes,
                         const bool print_warnings,
                         const std::string& constraint_filename = "",
                         const HyperedgeWeight constraint_weight = -100,
                         const HyperedgeWeight normal_edge_weight_bonus = 0,
                         OriginalEdgeSnapshot* out_snapshot = nullptr);

// Counts the (weighted) cut on the original hyperedge set captured in a
// snapshot. part_id must map HypernodeID -> PartitionID for the final
// partition.
HyperedgeWeight countOriginalCutHyperedges(const HyperedgeVector& original_hyperedges,
                                           const vec<HyperedgeWeight>& original_weights,
                                           const vec<PartitionID>& part_id);

// Counts the (weighted) cut on the original edge set (graph case) captured
// in a snapshot. part_id must map HypernodeID -> PartitionID for the final
// partition.
HyperedgeWeight countOriginalCutEdges(const EdgeVector& original_edges,
                                     const vec<HyperedgeWeight>& original_weights,
                                     const vec<PartitionID>& part_id);

// Convenience dispatcher: picks the graph or hypergraph cut counter based on
// snapshot.is_graph. Throws InvalidInputException if the snapshot was never
// populated (snapshot.valid == false).
HyperedgeWeight countOriginalCut(const OriginalEdgeSnapshot& snapshot,
                                 const vec<PartitionID>& part_id);

void addFixedVertices(mt_kahypar_hypergraph_t hypergraph,
                      const mt_kahypar_partition_id_t* fixed_vertices,
                      const PartitionID k);
void addFixedVerticesFromFile(mt_kahypar_hypergraph_t hypergraph,
                              const std::string& filename,
                              const PartitionID k);
void removeFixedVertices(mt_kahypar_hypergraph_t hypergraph);

}  // namespace io
}  // namespace mt_kahypar