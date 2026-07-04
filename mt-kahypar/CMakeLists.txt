add_subdirectory(partition)
add_subdirectory(parallel)
add_subdirectory(datastructures)
add_subdirectory(io)
add_subdirectory(utils)

target_precompile_headers(MtKaHyPar-Sources INTERFACE definitions.h)
target_precompile_headers(MtKaHyPar-Sources INTERFACE partition/context.h)

target_sources(MtKaHyPar-Sources INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/datastructures/dynamic_adjacency_array.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/datastructures/contraction_tree.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/datastructures/dynamic_graph_factory.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/datastructures/dynamic_graph.cpp
)
