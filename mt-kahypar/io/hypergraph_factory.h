#pragma once

#include <string>
#include <fstream> 

#include "include/mtkahypartypes.h"

#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/context_enum_classes.h"
#include "mt-kahypar/utils/cast.h"

namespace mt_kahypar {

// Forward declaration
class Context;

namespace io {

mt_kahypar_hypergraph_t readInputFile(const std::string& filename,
                                      const PresetType& preset,
                                      const InstanceType& instance,
                                      const FileFormat& format,
                                      const bool stable_construction,
                                      const bool remove_single_pin_hes,
                                      const bool print_warnings,
                                      const Context& context=Context(),
                                      mt_kahypar_hypergraph_t* original_snapshot = nullptr);

template<typename Hypergraph>
Hypergraph readInputFile(const std::string& filename,
                         const FileFormat& format,
                         const bool stable_construction,
                         const bool remove_single_pin_hes,
                         const bool print_warnings,
                         const Context& context = Context(),
                         mt_kahypar_hypergraph_t* original_snapshot = nullptr);

void addFixedVertices(mt_kahypar_hypergraph_t hypergraph,
                      const mt_kahypar_partition_id_t* fixed_vertices,
                      const PartitionID k);

void addFixedVerticesFromFile(mt_kahypar_hypergraph_t hypergraph,
                              const std::string& filename,
                              const PartitionID k);

void removeFixedVertices(mt_kahypar_hypergraph_t hypergraph);

}  // namespace io
}  // namespace mt_kahypar