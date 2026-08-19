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
  // multiplier avoids losing resolution to integer rounding once we
  // multiply by the fractional tuning_parameter below
  constexpr HyperedgeWeight multiplier = 10;
  constexpr HyperedgeWeight multiplier2 = 100;

  for (const HyperedgeID he : _hg.edges()) {
    const HyperedgeWeight original_weight = _hg.edgeWeight(he);
    const HyperedgeWeight original_scaled = original_weight * multiplier;
    _hhg.setEdgeWeight(he, original_scaled);

    if (original_weight < 0) {
      double min_positive_sum = std::numeric_limits<double>::max();
      bool found_any_pin_sum = false;

      for (const HypernodeID pin : _hg.pins(he)) {
        double accumulator = 0.0;
        for (const HyperedgeID incident_he : _hg.incidentEdges(pin)) {
          if (incident_he != he && _hg.edgeWeight(incident_he) > 0) {
            double partial_sum = static_cast<double>(_hg.edgeWeight(incident_he)) * multiplier;
            //if (_context.heuristicEdgeSize) {
            //  partial_sum *= (static_cast<double>(multiplier2) / _hg.edgeSize(incident_he));
            //}
            accumulator += partial_sum;
          }
        }
        if (accumulator < min_positive_sum) {
          min_positive_sum = accumulator;
          found_any_pin_sum = true;
        }
      }

      if (found_any_pin_sum) {
        const double heuristic_value = -min_positive_sum;
        const double t = _context.tuning_parameter; // expected in [0, 1]
	//LOG << "heuristic_value: " << heuristic_value;

        // Smooth interpolation:
        //   t = 0  -> original (scaled) weight
        //   t = 1  -> fully replaced by the heuristic value
        const double interpolated = (1.0 - t) * original_scaled + t * heuristic_value;
	//LOG << "iterpolated: " << interpolated;
  if (t < 0) {
    _hhg.setEdgeWeight(he, static_cast<HyperedgeWeight>(std::llround(t)));
  }
	else if (std::abs(interpolated) < std::abs(original_scaled)) {
        	_hhg.setEdgeWeight(he, static_cast<HyperedgeWeight>(std::llround(interpolated)));
	}
      }
      // if no pin has any positive incident edge, leave the original scaled weight in place
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
