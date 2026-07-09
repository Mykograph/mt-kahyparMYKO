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
    // Baseline: every edge is scaled by `multiplier` exactly once. This is the
    // one true "current scale" for _hhg — negative edges below either keep
    // this value or replace it, but never get multiplier applied a second time.
    const double scaled_original = _hg.edgeWeight(he) * multiplier;

    if (_hg.edgeWeight(he) >= 0) {
      _hhg.setEdgeWeight(he, static_cast<HyperedgeWeight>(std::llround(scaled_original)));
      continue;
    }

    // Negative edge: find the pin with the least positive support from its
    // *other* incident edges, and use that as an upper bound (in magnitude)
    // for how negative this edge is allowed to be.
    double weakest_pin_support = std::numeric_limits<double>::max();
    bool found_any_pin_with_support = false;

    for (const HypernodeID pin : _hg.pins(he)) {
      double accumulator = 0.0;
      bool pin_has_support = false;

      for (const HyperedgeID incident_he : _hg.incidentEdges(pin)) {
        if (incident_he != he && _hg.edgeWeight(incident_he) > 0) {
          double partial_sum = _hg.edgeWeight(incident_he) * multiplier;
          if (_context.heuristicEdgeSize) {
            // Real float ratio now, not integer division truncated to 0/1/2/...
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

    // If literally no pin has any positive support, there's nothing to bound
    // the negative edge by — fall back to leaving it at its original (scaled)
    // value rather than silently zeroing it out.
    if (!found_any_pin_with_support) {
      _hhg.setEdgeWeight(he, static_cast<HyperedgeWeight>(std::llround(scaled_original)));
      continue;
    }

    const double candidate = -weakest_pin_support * _context.tuning_parameter;

    // Both sides are now on the exact same scale (one factor of `multiplier`),
    // so this comparison is meaningful: "is the heuristic-bounded value less
    // negative (i.e. bigger) than the original scaled weight?"
    if (candidate > scaled_original) {
      _hhg.setEdgeWeight(he, static_cast<HyperedgeWeight>(std::llround(candidate)));
    } else {
      _hhg.setEdgeWeight(he, static_cast<HyperedgeWeight>(std::llround(scaled_original)));
    }
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
