#include <CLI/CLI.hpp>

#include <optional>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <random>
#include <unordered_set>

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"

#include "mt-kahypar/datastructures/hypergraph_common.h"

using Hypergraph = mt_kahypar::ds::StaticHypergraph;

using namespace mt_kahypar;

using namespace mt_kahypar;
namespace fs = std::filesystem;

void write_hgr_file(const fs::path output_path,
                    const io::HyperedgeVector& hyperedges,
                    const vec<HyperedgeWeight>& edge_weights,
                    const HypernodeID num_nodes) {
    std::ofstream out(output_path);
    
    HyperedgeID num_edges = hyperedges.size();
    
    // Header: weight_type 01 = nur edge weights
    out << num_edges << " " << num_nodes << " 01" << std::endl;
    
    // Hyperkanten mit Gewicht an erster Stelle
    for (HyperedgeID i = 0; i < num_edges; i++) {
        out << edge_weights[i] << " ";
        for (const HypernodeID pin : hyperedges[i]) {
            out << (pin + 1) << " ";  // +1 weil hMetis bei 1 anfängt
        }
        out << std::endl;
    }
    
    out.close();
}

void initialise_random_edge_weights(vec<HyperedgeWeight>& edge_weights, float p) {
    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<float> prob(0.0f, 1.0f);
    std::uniform_int_distribution<HyperedgeWeight> weight(2, 10);
    
    for (HyperedgeID i = 0; i < edge_weights.size(); i++) {
        if (prob(rng) < p) {
            edge_weights[i] = 1;//-weight(rng);  // random weight between 2 and 10
        } else {
            edge_weights[i] = 1;
        }
    }
}

HypernodeID generate_weights_from_hg(const fs::path hg_path,
                                      const fs::path weighted_path){

    HyperedgeID num_hyperedges = 0;
    HypernodeID num_hypernodes = 0;
    HyperedgeID num_removed = 0;
    io::HyperedgeVector hyperedges;
    vec<HyperedgeWeight> hyperedge_weights;
    vec<HypernodeWeight> hypernode_weights;

    io::readHypergraphFile(
        hg_path.string(),
        num_hyperedges,
        num_hypernodes,
        num_removed,
        hyperedges,
        hyperedge_weights,
        hypernode_weights,
        false,  // remove_single_pin_hes
        true    // print_warnings
    );

    


 

    bool has_edge_weights = !hyperedge_weights.empty();
    bool has_node_weights = !hypernode_weights.empty();
    std::cout << "Has edge weights: " << std::to_string(has_edge_weights) << std::endl;
    std::cout << "Has node weights: " << std::to_string(has_node_weights) << std::endl;
    if (has_edge_weights || has_node_weights) {
        std::cout << "already has edge weights or node weights" << std::endl;
        return num_hypernodes;
    }

    //initialise the edge weights with 1's
    vec<HyperedgeWeight> edge_weights(num_hyperedges);
    initialise_random_edge_weights(edge_weights, 0.3f); // 30% chance to have a weight > 1

    /*for (HyperedgeID i = 0; i < num_hyperedges; i++) {
        edge_weights[i] = 1;//pick_random(10) + 1; // random weight between 1 and 10
    }*/

    //write to file
    write_hgr_file(weighted_path, hyperedges, edge_weights, num_hypernodes);
    return num_hypernodes;
}

std::vector<mt_kahypar::PartitionID> readPartitionFile(const fs::path partition_path, const HypernodeID num_nodes) {
    std::vector<mt_kahypar::PartitionID> partition(num_nodes);
    std::ifstream in(partition_path);
    if (!in.is_open()) {
        throw std::runtime_error("Could not open partition file: " + partition_path.string());
    }
    for (HypernodeID i = 0; i < num_nodes; i++) {
        if (!(in >> partition[i])) {
            throw std::runtime_error("Error reading partition file at line " + std::to_string(i + 1));
        }
    }
    in.close();
    return partition;
}


void setEdgeWeightsBasedOnPartition(
    const io::HyperedgeVector& hyperedges,
    const std::vector<PartitionID>& partition,
    vec<HyperedgeWeight>& edge_weights)
{
    for (HyperedgeID i = 0; i < hyperedges.size(); i++) {
        const auto& he = hyperedges[i];
        if (he.empty()) { 
            std::cout << "Warning: Hyperedge " << i << " is empty. Assigning default weight of 1." << std::endl;
            edge_weights[i] = 1; continue; }
        PartitionID first = partition[he[0]];
        bool is_cut_edge = false;
        for (HypernodeID pin : he) {
            if (partition[pin] != first) { is_cut_edge = true; break; }
        }
        edge_weights[i] = is_cut_edge ? -10 : 1;
        if(is_cut_edge) {
            std::cout << "Hyperedge " << i << " is a cut edge. Assigning weight -10." << std::endl;
        }
    }
}

HypernodeID generate_weights_from_hgp(const fs::path hg_path,
                                       const fs::path weighted_path,
                                       const fs::path partition_path,
                                       int k)
{
    HyperedgeID num_hyperedges = 0;
    HypernodeID num_hypernodes = 0;
    HyperedgeID num_removed    = 0;
    io::HyperedgeVector hyperedges;
    vec<HyperedgeWeight> hyperedge_weights;
    vec<HypernodeWeight> hypernode_weights;

    io::readHypergraphFile(
        hg_path.string(),
        num_hyperedges, num_hypernodes, num_removed,
        hyperedges, hyperedge_weights, hypernode_weights,
        false, true
    );

    // Weight-Check VOR allem anderen
    if (!hyperedge_weights.empty() || !hypernode_weights.empty()) {
        std::cout << "already has edge weights or node weights, will be replaced" << std::endl;
    }

    // Partition lesen
    std::vector<PartitionID> partition = readPartitionFile(partition_path, num_hypernodes);

    // Kantengewichte direkt aus io::HyperedgeVector berechnen - kein HG-Objekt nötig
    vec<HyperedgeWeight> edge_weights(num_hyperedges);
    setEdgeWeightsBasedOnPartition(hyperedges, partition, edge_weights);

    write_hgr_file(weighted_path, hyperedges, edge_weights, num_hypernodes);
    return num_hypernodes;
}


int main(int argc, char* argv[]) {
    bool partitioned_provided = false;
    std::string hypergraph_path;
    std::string weighted_dir;

    std::string partition_path;
    int k = 2;

    CLI::App app;

    //required arguments
    app.add_option("-H,--hypergraph", hypergraph_path, "Hypergraph Filename or directory")->required();
    app.add_option("-w,--weighted", weighted_dir, "Weighted directory")->required();

    //optional arguments
    CLI::Option* opt = app.add_option("-p,--partitioned", partition_path, "Partition path (optional)");
    CLI::Option* opt_k = app.add_option("-k,--k", k, "Number of partitions (optional) - only needed if partition is provided)")->needs(opt);

    CLI11_PARSE(app, argc, argv);
     

    if (opt->count() > 0) {
        partitioned_provided = true;
    }

    // dann statt fs::path einfach:
    fs::path weighted_dir_path = weighted_dir;
    fs::path hypergraph_path_p = hypergraph_path;

    std::cout << "Generating weights for hypergraph: " << hypergraph_path_p << std::endl;

    if (!fs::is_directory(weighted_dir_path)) {
            throw fs::filesystem_error(
                "weighted path is not a directory",
                weighted_dir_path,
                std::make_error_code(std::errc::not_a_directory)
            );
        }

    if (!(fs::is_regular_file(hypergraph_path_p) && hypergraph_path_p.extension() == ".hgr")) {
        throw fs::filesystem_error(
            "hypergraph path is not a .hgr file",
            hypergraph_path_p,
            std::make_error_code(std::errc::invalid_argument)
        );
    }

    if (!(fs::is_regular_file(hypergraph_path_p))) {
        throw fs::filesystem_error(
            "partition path is not a file",
            partition_path,
            std::make_error_code(std::errc::invalid_argument)
        );
    }

    
    HypernodeID generated_weights;
    
    if (opt->count() > 0) {
        fs::path weighted_file = weighted_dir_path / (hypergraph_path_p.filename().string() + "-withWeights-withPartition" + ".hgr");
        generated_weights = generate_weights_from_hgp(hypergraph_path_p, weighted_file, partition_path, k);
    } else {
        fs::path weighted_file = weighted_dir_path / (hypergraph_path_p.filename().string() + "-withWeights" + ".hgr");
        generated_weights = generate_weights_from_hg(hypergraph_path_p, weighted_file);
    }
    return 0;
}
