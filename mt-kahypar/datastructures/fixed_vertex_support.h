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

#include <atomic>
#include <memory>
#include <unordered_map>

#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"

namespace mt_kahypar {
namespace ds {

// Forward declaration to avoid circular include with dynamic_graph.h
class DynamicGraph;

template<class Hypergraph>
class FixedVertexSupport {

  static constexpr bool debug = false;

  struct FixedVertexData {
    // ! Fixed vertex block ID
    PartitionID block;
    // ! Number of fixed vertices contracted onto this node
    HypernodeID fixed_vertex_contraction_cnt;
    // ! Weight at the time it becomes fixed
    HypernodeWeight fixed_vertex_weight;
    // ! Spin lock to syncronize contractions
    SpinLock sync;
  };

 public:
  FixedVertexSupport();

  FixedVertexSupport(const HypernodeID num_nodes,
                     const PartitionID k);

  FixedVertexSupport(const FixedVertexSupport&) = delete;
  FixedVertexSupport& operator=(const FixedVertexSupport&) = delete;

  // Declared but defined in .cpp to allow unique_ptr<DynamicGraph> with forward declaration
  FixedVertexSupport(FixedVertexSupport&&);
  FixedVertexSupport& operator=(FixedVertexSupport&&);
  ~FixedVertexSupport();

  void setHypergraph(const Hypergraph* hg) {
    _hg = hg;
  }

  void setMaxBlockWeight(const std::vector<HypernodeWeight> max_block_weights) {
    if ( hasFixedVertices() ) {
      ASSERT(max_block_weights.size() >= static_cast<size_t>(_k));
      _max_block_weights = max_block_weights;
    }
  }

  PartitionID numBlocks() const {
    return _k;
  }

  // ####################### Fixed Vertex Block Weights #######################

  bool hasFixedVertices() const {
    return _total_fixed_vertex_weight.load(std::memory_order_relaxed) > 0;
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HypernodeWeight totalFixedVertexWeight() const {
    return _total_fixed_vertex_weight.load(std::memory_order_relaxed);
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE HypernodeWeight fixedVertexBlockWeight(const PartitionID block) const {
    ASSERT(block != kInvalidPartition && block < _k);
    return _fixed_vertex_block_weights[block].load(std::memory_order_relaxed);
  }

  // ####################### Fixed Vertex Information #######################

  void fixToBlock(const HypernodeID hn, const PartitionID block) {
    ASSERT(_hg);
    ASSERT(hn < _num_nodes);
    ASSERT(block != kInvalidPartition && block < _k);
    PartitionID expected = kInvalidPartition;
    PartitionID desired = block;
    if ( __atomic_compare_exchange_n(&_fixed_vertex_data[hn].block,
           &expected, desired, false, __ATOMIC_ACQ_REL, __ATOMIC_RELAXED) ) {
      const HypernodeWeight weight_of_hn = _hg->nodeWeight(hn);
      _fixed_vertex_data[hn].fixed_vertex_contraction_cnt = 1;
      _fixed_vertex_data[hn].fixed_vertex_weight = weight_of_hn;
      _fixed_vertex_block_weights[block].fetch_add(
        weight_of_hn, std::memory_order_relaxed);
      _total_fixed_vertex_weight.fetch_add(
        weight_of_hn, std::memory_order_relaxed);
    } else {
      ASSERT(_fixed_vertex_data[hn].block == block,
        "Try to fix hypernode" << hn << "to block" << block
        << ", but it is already fixed to block" << _fixed_vertex_data[hn].block);
    }
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE bool isFixed(const HypernodeID hn) const {
    return hn < _num_nodes && fixedVertexBlock(hn) != kInvalidPartition;
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE PartitionID fixedVertexBlock(const HypernodeID hn) const {
    ASSERT(hn < _num_nodes);
    return __atomic_load_n(&_fixed_vertex_data[hn].block, __ATOMIC_RELAXED);
  }

  // ####################### (Un)contractions #######################

  bool contract(const HypernodeID u, const HypernodeID v);
  bool contractWithoutChains(const HypernodeID u, const HypernodeID v);
  void uncontract(const HypernodeID u, const HypernodeID v);

  // ####################### Negative Constraints #######################

  bool hasNegativeConstraints() const {
    return _constraint_graph != nullptr;
  }

  // ! Builds the constraint graph from the given pairs of node IDs.
  // ! Defined in fixed_vertex_support.cpp to avoid circular includes.
  void setNegativeConstraints(const vec<std::pair<HypernodeID, HypernodeID>>& constraints);

  // ! Returns the constraint graph. Only valid if hasNegativeConstraints() is true.
  const DynamicGraph& getConstraintGraph() const {
    ASSERT(_constraint_graph != nullptr);
    return *_constraint_graph;
  }

  // ! Maps a hypergraph node ID to its ID in the constraint graph.
  // ! Returns false if the node is not involved in any constraint.
  bool getConstraintIdFromHypergraphId(const HypernodeID hg_id, HypernodeID& out_id) const {
    auto it = _hg_id_to_constraint_id.find(hg_id);
    if (it == _hg_id_to_constraint_id.end()) return false;
    out_id = it->second;
    return true;
  }

  // ####################### Miscellaneous #######################

  FixedVertexSupport<Hypergraph> copy() const {
    FixedVertexSupport<Hypergraph> cpy;
    cpy._num_nodes = _num_nodes;
    cpy._k = _k;
    cpy._hg = _hg;
    cpy._total_fixed_vertex_weight = _total_fixed_vertex_weight;
    cpy._fixed_vertex_block_weights = _fixed_vertex_block_weights;
    cpy._max_block_weights = _max_block_weights;
    cpy._fixed_vertex_data = _fixed_vertex_data;
    // _constraint_graph and _hg_id_to_constraint_id not copied:
    // they are read-only after construction and large; caller can rebuild if needed
    return cpy;
  }

  size_t size_in_bytes() const {
    return (sizeof(CAtomic<HypernodeWeight>) + sizeof(HypernodeWeight)) * _k +
      sizeof(FixedVertexData) * _num_nodes;
  }

 private:
  bool contractImpl(const HypernodeID u, const HypernodeID v, bool ignore_v);

  // ! Number of nodes
  HypernodeID _num_nodes;
  // ! Number of blocks
  PartitionID _k;
  // ! Underlying hypergraph
  const Hypergraph* _hg;
  // ! Total weight of all fixed vertices
  CAtomic<HypernodeWeight> _total_fixed_vertex_weight;
  // ! Weight of all vertices fixed to a block
  vec<CAtomic<HypernodeWeight>> _fixed_vertex_block_weights;
  // ! Maximum allowed fixed vertex block weight
  std::vector<HypernodeWeight> _max_block_weights;
  // ! Fixed vertex block IDs of each node
  vec<FixedVertexData> _fixed_vertex_data;
  // ! Graph representing the negative constraints (node weight = original hypergraph ID)
  std::unique_ptr<DynamicGraph> _constraint_graph;
  // ! Map from hypergraph node ID to constraint graph node ID
  std::unordered_map<HypernodeID, HypernodeID> _hg_id_to_constraint_id;
};

}  // namespace ds
}  // namespace mt_kahypar
