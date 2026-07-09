/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/refinement/i_refiner.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/partition/coarsening/coarsening_commons.h"
#include "mt-kahypar/utils/timer.h"


namespace mt_kahypar {

template<typename TypeTraits>
class MultilevelCoarsenerBase {
 private:

  static constexpr bool debug = false;
  using Hypergraph = typename TypeTraits::Hypergraph;
  using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;

 public:
  MultilevelCoarsenerBase(Hypergraph& hypergraph,
                          const Context& context,
                          UncoarseningData<TypeTraits>& uncoarseningData) :
          _hg(hypergraph),
          _context(context),
          _timer(utils::Utilities::instance().getTimer(context.utility_id)),
          _uncoarseningData(uncoarseningData),
          _hhg(hypergraph.copy()) {
            //_hhg.setEdgeWeight(1, 999);
            //LOG << "_hg.edgeWeight(1): " << _hg.edgeWeight(1);
            //LOG << "_hhg.edgeWeight(1): " << _hhg.edgeWeight(1);
          heuristicHypergraph(); 
          _uncoarseningData.setHeuristicHypergraph(_hhg);    
}

  MultilevelCoarsenerBase(const MultilevelCoarsenerBase&) = delete;
  MultilevelCoarsenerBase(MultilevelCoarsenerBase&&) = delete;
  MultilevelCoarsenerBase & operator= (const MultilevelCoarsenerBase &) = delete;
  MultilevelCoarsenerBase & operator= (MultilevelCoarsenerBase &&) = delete;

  virtual ~MultilevelCoarsenerBase() = default;

 protected:

  HypernodeID currentNumNodes() const {
    if ( _uncoarseningData.hierarchy.empty() ) {
      return _hg.initialNumNodes();
    } else {
      return _uncoarseningData.hierarchy.back().contractedHypergraph().initialNumNodes();
    }
  }

  Hypergraph& currentHypergraph() {
    if ( _uncoarseningData.hierarchy.empty() ) {
      //LOG << "Returning original hypergraph *************************************************";
      //LOG << "original num nodes: " << _hg.initialNumNodes();
      return _hhg;
    } else {
      //LOG << "Returning coarsened hypergraph *************************************************";
      return _uncoarseningData.hierarchy.back().contractedHypergraph();
    }
  }

void heuristicHypergraph() {
  // multiplier keeps enough precision when later multiplying by tuning_parameter (in [0,1])
  // before we round back down to an integral HyperedgeWeight.
  constexpr double multiplier = 10.0;
  constexpr double multiplier2 = 100.0;

  for (const HyperedgeID he : _hg.edges()) {
    // Baseline: every edge is scaled by `multiplier` exactly once.
    const double scaled_original = _hg.edgeWeight(he) * multiplier;

    if (_hg.edgeWeight(he) >= 0) {
      _hhg.setEdgeWeight(he, static_cast<HyperedgeWeight>(std::llround(scaled_original)));
      continue;
    }

    // Negative edge: for each pin, sum the weights of its *other* positive
    // incident edges. The heuristic bound is the minimum such sum over all
    // pins of he — i.e. the "worst supported" pin caps how negative he is
    // allowed to be.
    double weakest_pin_support = std::numeric_limits<double>::max();
    bool found_any_pin_with_support = false;

    for (const HypernodeID pin : _hg.pins(he)) {
      double accumulator = 0.0;
      bool pin_has_support = false;

      for (const HyperedgeID incident_he : _hg.incidentEdges(pin)) {
        if (incident_he != he && _hg.edgeWeight(incident_he) > 0) {
          double partial_sum = _hg.edgeWeight(incident_he) * multiplier;
          if (_context.heuristicEdgeSize) {
            partial_sum *= (multiplier2 / static_cast<double>(_hg.edgeSize(incident_he)));
          }
          accumulator += partial_sum;
          pin_has_support = true;
        }
      }

      if (pin_has_support && accumulator < weakest_pin_support) {
        weakest_pin_support = accumulator;
        found_any_pin_with_support = true;
      }
    }

    // No pin has any positive support at all -> nothing to blend towards,
    // keep the edge as-is.
    if (!found_any_pin_with_support) {
      _hhg.setEdgeWeight(he, static_cast<HyperedgeWeight>(std::llround(scaled_original)));
      continue;
    }

    // Heuristic bound: negative, on the same scale as scaled_original.
    const double heuristic_bound = -weakest_pin_support;

    // Smooth blend: t=0 -> untouched original, t=1 -> full heuristic bound.
    // Monotonic and continuous in tuning_parameter, no hard switch.
    const double t = _context.tuning_parameter;
    const double blended = (1.0 - t) * scaled_original + t * heuristic_bound;

    _hhg.setEdgeWeight(he, static_cast<HyperedgeWeight>(std::llround(blended)));
  }
}

  PartitionedHypergraph& currentPartitionedHypergraph() {
    ASSERT(_uncoarseningData.is_finalized);
    return *_uncoarseningData.partitioned_hg;
   }

 protected:
  Hypergraph& _hg;
  Hypergraph _hhg;
  const Context& _context;
  utils::Timer& _timer;
  UncoarseningData<TypeTraits>& _uncoarseningData;
};
}  // namespace mt_kahypar
