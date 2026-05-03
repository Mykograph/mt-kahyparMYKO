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
namespace fs = std::filesystem;

const int DEFAULT_GENERATE_TYPE = 0; // 0: random weights, 1: partition-based weights
const float DEFAULT_RANDOM_WEIGHT_PROBABILITY = 0.3f; // 30% chance to have a weight > 1
const int NEGATIVE_WEIGHT = 10; // Penalty for cut edges in partition-based weighting
const int INVERSE_WEIGHT_MULTIPLIER = 10; // Multiplier for inverse weights in partition-based weighting

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
            edge_weights[i] = -weight(rng);  // random weight between 2 and 10
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
    std::string* output,
    int mode)
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
        if(mode == 0) {
            edge_weights[i] = is_cut_edge ? -NEGATIVE_WEIGHT : edge_weights[i]; // cut edges get weight -10, others keep their weight (or get 1 if not previously set)
        } else if (mode == 1) {
            edge_weights[i] = is_cut_edge ? edge_weights[i] - NEGATIVE_WEIGHT : edge_weights[i]; // cut edges get weight 10, others keep their weight (or get 1 if not previously set)
        } else {
            edge_weights[i] = is_cut_edge ? -(int)( (double)1 / (double) edge_weights[i] * INVERSE_WEIGHT_MULTIPLIER) : edge_weights[i]; // cut edges get inverse weight, others keep their weight (or get 1 if not previously set)
        }
        
        if(is_cut_edge) {
            *output += ("Hyperedge " + std::to_string(i) + " is a cut edge. Assigning weight -10.\n");
        }
    }
}

HypernodeID generate_weights_from_hgp(const fs::path hg_path,
                                       const fs::path weighted_path,
                                       const fs::path partition_path,
                                       int generate_type,
                                       bool output_print)
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

    if (partition_path.empty()) {
        switch (generate_type) {
            case 0:
                std::fill(edge_weights.begin(), edge_weights.end(), 1);
                break;
            case 1:
        
                initialise_random_edge_weights(edge_weights, DEFAULT_RANDOM_WEIGHT_PROBABILITY ); 
        }
    } else {
        // Partition lesen
        std::vector<PartitionID> partition = readPartitionFile(partition_path, num_hypernodes);
        // Falls bereits edge weights vorhanden sind, diese als Basis nehmen, ansonsten mit 1 initialisieren
                if (hyperedge_weights.empty()) {
                    std::fill(edge_weights.begin(), edge_weights.end(), 1);
                } else {
                    edge_weights = hyperedge_weights;
                }  
        switch (generate_type) {
            case 0:
                setEdgeWeightsBasedOnPartition(hyperedges, partition, edge_weights, &output, 0);
                break;
            case 1:
                setEdgeWeightsBasedOnPartition(hyperedges, partition, edge_weights, &output, 1);
                break;
            case 2: 
                setEdgeWeightsBasedOnPartition(hyperedges, partition, edge_weights, &output, 2);
        }
    }   
    write_hgr_file(weighted_path, hyperedges, edge_weights, num_hypernodes);

    if (!output_print) {
        std::cout << output;
    }
    return num_hypernodes;
}


int main(int argc, char* argv[]) {
    std::string hypergraph_path;
    std::string weighted_dir;

    std::string partition_path;
    int generate_type = 0; 

    CLI::App app;

    //required arguments
    app.add_option("-H,--hypergraph", hypergraph_path, "Hypergraph Filename or directory")->required();
    app.add_option("-w,--weighted", weighted_dir, "Weighted directory")->required();

    app.add_option("-g,--generate-type", generate_type, "Type of weight generation: 0 for random, 1 for partition-based (optional)")->check(CLI::Range(0, 1));

    //optional arguments
    CLI::Option* opt = app.add_option("-p,--partitioned", partition_path, "Partition path (optional)");

    bool output_print = false;
    app.add_flag("-n,--no-output", output_print, "Disable output to console (optional)")->default_val("false")->capture_default_str();

        

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

    if (!(fs::is_regular_file(hypergraph_path_p))) {
        throw fs::filesystem_error(
            "partition path is not a file",
            partition_path,
            std::make_error_code(std::errc::invalid_argument)
        );
    }

    
    HypernodeID generated_weights;
    
    // Generate weights based on the provided options
    //Optiona without partition
    //option 0: random weights,
    //option 1: all 1's
    //
    //Optiona with partition
    //option 0: weights based on partition
    //option 1: weights based on partition with relative weights (cut edges weight = weight -x)
    //option 2: weights based on partition with relative weights (cut edges weight = weight * -(1/x))
    

    if ((fs::is_regular_file(hypergraph_path_p) && hypergraph_path_p.extension() == ".hgr")) {
        if (opt->count() > 0) {
        fs::path weighted_file = weighted_dir_path / (hypergraph_path_p.filename().string() + "-withWeights-withPartition:"+ partition_path_p.filename().string() +  ".hgr");
        generated_weights = generate_weights_from_hgp(hypergraph_path_p, weighted_file, partition_path_p, generate_type, output_print);
    } else {
        fs::path weighted_file = weighted_dir_path / (hypergraph_path_p.filename().string() + "-withWeights" + ".hgr");
        generated_weights = generate_weights_from_hgp(hypergraph_path_p, weighted_file, fs::path(), generate_type, output_print);
    }
    } else if (fs::is_directory(hypergraph_path_p)) {
        for (const auto& entry : fs::directory_iterator(hypergraph_path_p)) {
            if (entry.is_regular_file() && entry.path().extension() == ".hgr") {
                fs::path weighted_file;
                if (opt->count() > 0) {
                    weighted_file = weighted_dir_path / (entry.path().filename().string() + "-withWeights-withPartition:"+ partition_path_p.filename().string() +  ".hgr");
                    generated_weights = generate_weights_from_hgp(entry.path(), weighted_file, partition_path_p, generate_type, output_print);
                } else {
                    weighted_file = weighted_dir_path / (entry.path().filename().string() + "-withWeights" + ".hgr");
                    generated_weights = generate_weights_from_hgp(entry.path(), weighted_file, fs::path(), generate_type, output_print);
                }
            }
        }
    } else {
        std::cerr << "Invalid hypergraph path: " << hypergraph_path_p << std::endl;
        return 1;
    }

    
    return 0;
}
