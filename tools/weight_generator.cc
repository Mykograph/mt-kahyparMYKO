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

void output_function(const std::string& message) {
    std::cout << "#############Finished Program################" << std::endl;
    std::cout << message << std::endl;
}

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
    vec<HyperedgeWeight>& edge_weights,
    std::string* output)
{
    for (HyperedgeID i = 0; i < hyperedges.size(); i++) {
        const auto& he = hyperedges[i];
        if (he.empty()) { 
            *output += ("Warning: Hyperedge " + std::to_string(i) + " is empty. Assigning default weight of 1.\n");
            edge_weights[i] = 1; continue; }
        PartitionID first = partition[he[0]];
        bool is_cut_edge = false;
        for (HypernodeID pin : he) {
            if (partition[pin] != first) { is_cut_edge = true; break; }
        }
        edge_weights[i] = is_cut_edge ? -10 : 1;
        if(is_cut_edge) {
            *output += ("Hyperedge " + std::to_string(i) + " is a cut edge. Assigning weight -10.\n");
        }
    }
}

HypernodeID generate_weights_from_hgp(const fs::path hg_path,
                                       const fs::path weighted_path,
                                       const fs::path partition_path,
                                       int k=2,
                                       int generate_type = 0)
{
    HyperedgeID num_hyperedges = 0;
    HypernodeID num_hypernodes = 0;
    HyperedgeID num_removed    = 0;
    io::HyperedgeVector hyperedges;
    vec<HyperedgeWeight> hyperedge_weights;
    vec<HypernodeWeight> hypernode_weights;
    std::string output = "";

    io::readHypergraphFile(
        hg_path.string(),
        num_hyperedges, num_hypernodes, num_removed,
        hyperedges, hyperedge_weights, hypernode_weights,
        false, true
    );

    if (!hyperedge_weights.empty() || !hypernode_weights.empty()) {
        output += "Warning: The input hypergraph file " + hg_path.string() + " already has edge weights or node weights, these will be replaced. \n";
    }
    
            vec<HyperedgeWeight> edge_weights(num_hyperedges);

    switch (generate_type) {
        case 0:
            
            initialise_random_edge_weights(edge_weights, 0.3f ); // 30% chance to have a weight > 1

            write_hgr_file(weighted_path, hyperedges, edge_weights, num_hypernodes);
            return num_hypernodes;
        case 1:
            // Partition lesen
            std::vector<PartitionID> partition = readPartitionFile(partition_path, num_hypernodes);
            setEdgeWeightsBasedOnPartition(hyperedges, partition, edge_weights, &output);

            write_hgr_file(weighted_path, hyperedges, edge_weights, num_hypernodes);
            return num_hypernodes;
    }

    std::cout << output;
}


int main(int argc, char* argv[]) {
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

    fs::path weighted_dir_path = weighted_dir;
    fs::path hypergraph_path_p = hypergraph_path;
    fs::path partition_path_p = partition_path;


    // Check if paths are valid

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
        fs::path weighted_file = weighted_dir_path / (hypergraph_path_p.filename().string() + "-withWeights-withPartition:"+ partition_path_p.filename().string() + ".hgr");
        generated_weights = generate_weights_from_hgp(hypergraph_path_p, weighted_file, partition_path_p, k);
    } else {
        fs::path weighted_file = weighted_dir_path / (hypergraph_path_p.filename().string() + "-withWeights" + ".hgr");
        generated_weights = generate_weights_from_hgp(hypergraph_path_p, weighted_file, fs::path(), k, 0);
    }
    return 0;
}
