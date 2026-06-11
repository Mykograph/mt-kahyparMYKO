#include <CLI/CLI.hpp>

#include <optional>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <random>
#include <unordered_set>

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"

using Hypergraph = mt_kahypar::ds::StaticHypergraph;

using namespace mt_kahypar;
namespace fs = std::filesystem;

const int DEFAULT_GENERATE_TYPE = 0;
const float DEFAULT_RANDOM_WEIGHT_PROBABILITY = 0.3f;
int NEGATIVE_WEIGHT = 10;
const int INVERSE_WEIGHT_MULTIPLIER = 10;

// ============================================================
// Shared utility
// ============================================================

void output_function(const std::string& message) {
    std::cout << "#############Finished Program################" << std::endl;
    std::cout << message << std::endl;
}

std::vector<mt_kahypar::PartitionID> readPartitionFile(const fs::path& partition_path,
                                                        const HypernodeID num_nodes)
{
    std::vector<mt_kahypar::PartitionID> partition(num_nodes);
    std::ifstream in(partition_path);
    if (!in.is_open())
        throw std::runtime_error("Could not open partition file: " + partition_path.string());
    for (HypernodeID i = 0; i < num_nodes; i++) {
        if (!(in >> partition[i]))
            throw std::runtime_error("Error reading partition file at line " + std::to_string(i + 1));
    }
    return partition;
}

// ============================================================
// Hypergraph functions (unchanged from original)
// ============================================================

void write_hgr_file(const fs::path output_path,
                    const io::HyperedgeVector& hyperedges,
                    const vec<HyperedgeWeight>& edge_weights,
                    const HypernodeID num_nodes)
{
    std::ofstream out(output_path);
    HyperedgeID num_edges = hyperedges.size();
    out << num_edges << " " << num_nodes << " 01" << std::endl;
    for (HyperedgeID i = 0; i < num_edges; i++) {
        out << edge_weights[i] << " ";
        for (const HypernodeID pin : hyperedges[i])
            out << (pin + 1) << " ";
        out << std::endl;
    }
}

void initialise_random_edge_weights(vec<HyperedgeWeight>& edge_weights, float p) {
    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<float> prob(0.0f, 1.0f);
    std::uniform_int_distribution<HyperedgeWeight> weight(2, 10);
    for (HyperedgeID i = 0; i < edge_weights.size(); i++) {
        if (prob(rng) < p)
            edge_weights[i] = -weight(rng);
        else
            edge_weights[i] = 1;
    }
}

void setEdgeWeightsBasedOnPartition(
    const io::HyperedgeVector& hyperedges,
    const std::vector<PartitionID>& partition,
    vec<HyperedgeWeight>& edge_weights,
    std::string* output,
    int mode,
    double negative_weight_chance)
{
    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<float> prob(0.0f, 1.0f);

    for (HyperedgeID i = 0; i < hyperedges.size(); i++) {
        const auto& he = hyperedges[i];
        if (he.empty()) {
            *output += "Warning: Hyperedge " + std::to_string(i) + " is empty. Assigning default weight of 1.\n";
            edge_weights[i] = 1;
            continue;
        }
        PartitionID first = partition[he[0]];
        bool is_cut_edge = false;
        for (HypernodeID pin : he) {
            if (partition[pin] != first) { is_cut_edge = true; break; }
        }
        if (mode == 0) {
            edge_weights[i] = is_cut_edge ? -NEGATIVE_WEIGHT : edge_weights[i];
        } else if (mode == 1) {
            edge_weights[i] = is_cut_edge ? edge_weights[i] - NEGATIVE_WEIGHT : edge_weights[i];
        } else if (mode == 2) {
            edge_weights[i] = is_cut_edge
                ? -(int)((double)1 / (double)edge_weights[i] * INVERSE_WEIGHT_MULTIPLIER)
                : edge_weights[i];
        } else {
            edge_weights[i] = is_cut_edge
                ? (prob(rng) < negative_weight_chance ? -NEGATIVE_WEIGHT : edge_weights[i])
                : edge_weights[i];
        }
        if (is_cut_edge)
            *output += "Hyperedge " + std::to_string(i) + " is a cut edge. Assigning weight -10.\n";
    }
}

HypernodeID generate_weights_from_hgp(const fs::path hg_path,
                                       const fs::path weighted_path,
                                       const fs::path partition_path,
                                       int generate_type,
                                       bool output_print,
                                       double stochastic)
{
    HyperedgeID num_hyperedges = 0;
    HypernodeID num_hypernodes = 0;
    HyperedgeID num_removed    = 0;
    io::HyperedgeVector hyperedges;
    vec<HyperedgeWeight> hyperedge_weights;
    vec<HypernodeWeight> hypernode_weights;
    std::string output;

    io::readHypergraphFile(
        hg_path.string(),
        num_hyperedges, num_hypernodes, num_removed,
        hyperedges, hyperedge_weights, hypernode_weights,
        false, true
    );

    if (!hyperedge_weights.empty() || !hypernode_weights.empty())
        output += "Warning: The input hypergraph file " + hg_path.string()
               + " already has edge weights or node weights, these will be replaced.\n";

    vec<HyperedgeWeight> edge_weights(num_hyperedges);

    if (partition_path.empty()) {
        switch (generate_type) {
            case 0: std::fill(edge_weights.begin(), edge_weights.end(), 1); break;
            case 1: initialise_random_edge_weights(edge_weights, DEFAULT_RANDOM_WEIGHT_PROBABILITY); break;
        }
    } else {
        std::vector<PartitionID> partition = readPartitionFile(partition_path, num_hypernodes);
        if (hyperedge_weights.empty())
            std::fill(edge_weights.begin(), edge_weights.end(), 1);
        else
            edge_weights = hyperedge_weights;

        switch (generate_type) {
            case 0: setEdgeWeightsBasedOnPartition(hyperedges, partition, edge_weights, &output, 0, 0.0); break;
            case 1: setEdgeWeightsBasedOnPartition(hyperedges, partition, edge_weights, &output, 1, 0.0); break;
            case 2: setEdgeWeightsBasedOnPartition(hyperedges, partition, edge_weights, &output, 2, 0.0); break;
            case 3: setEdgeWeightsBasedOnPartition(hyperedges, partition, edge_weights, &output, 3, stochastic); break;
        }
    }

    write_hgr_file(weighted_path, hyperedges, edge_weights, num_hypernodes);

    if (!output_print)
        std::cout << output;

    return num_hypernodes;
}

// ============================================================
// Plain graph (METIS .graph / .metis) support
// ============================================================

struct GraphEdge {
    int to;
    int weight;
};

using AdjacencyList = std::vector<std::vector<GraphEdge>>;

// Read a METIS graph file.
// Header: <num_vertices> <num_edges> [fmt] [ncon]
//   fmt ones digit  = 1 → edge weights present
//   fmt tens digit  = 1 → vertex weights present (skipped)
AdjacencyList read_metis_graph(const fs::path& graph_path,
                                int& num_nodes,
                                int& num_edges_total,
                                bool& had_edge_weights)
{
    std::ifstream in(graph_path);
    if (!in.is_open())
        throw std::runtime_error("Could not open graph file: " + graph_path.string());

    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line[0] != '%') break;
    }

    std::istringstream header(line);
    int fmt = 0, ncon = 1;
    header >> num_nodes >> num_edges_total >> fmt >> ncon;

    bool has_edge_weights   = (fmt % 10) == 1;
    bool has_vertex_weights = ((fmt / 10) % 10) == 1;
    had_edge_weights = has_edge_weights;

    AdjacencyList adj(num_nodes);
    int node_idx = 0;

    while (std::getline(in, line) && node_idx < num_nodes) {
        if (line.empty() || line[0] == '%') continue;
        std::istringstream ss(line);

        // Skip vertex weights (ncon values per node)
        if (has_vertex_weights) {
            for (int c = 0; c < ncon; ++c) {
                int vw; ss >> vw;
            }
        }

        if (has_edge_weights) {
            int neighbor, ew;
            while (ss >> neighbor >> ew)
                adj[node_idx].push_back({neighbor - 1, ew}); // → 0-indexed
        } else {
            int neighbor;
            while (ss >> neighbor)
                adj[node_idx].push_back({neighbor - 1, 1});
        }
        ++node_idx;
    }
    return adj;
}

// Write a weighted METIS graph (fmt=1, edge weights).
// edge_weights_map[u][i] is the weight for the i-th neighbor of node u.
void write_metis_graph(const fs::path& output_path,
                       const AdjacencyList& adj,
                       const std::vector<std::vector<int>>& edge_weights_map)
{
    long long total_half_edges = 0;
    for (const auto& neighbors : adj)
        total_half_edges += neighbors.size();
    long long num_edges = total_half_edges / 2;

    std::ofstream out(output_path);
    if (!out.is_open())
        throw std::runtime_error("Could not open output path: " + output_path.string());

    int num_nodes = (int)adj.size();
    out << num_nodes << " " << num_edges << " 1\n";

    for (int u = 0; u < num_nodes; ++u) {
        for (int i = 0; i < (int)adj[u].size(); ++i) {
            if (i > 0) out << " ";
            out << (adj[u][i].to + 1) << " " << edge_weights_map[u][i];
        }
        out << "\n";
    }
}

// Assign weights to graph edges, mirroring the hypergraph logic.
std::vector<std::vector<int>> assign_graph_edge_weights(
    const AdjacencyList& adj,
    const std::vector<PartitionID>& partition,   // empty = no partition
    int generate_type,
    double stochastic,
    std::string& output)
{
    int num_nodes = (int)adj.size();
    std::vector<std::vector<int>> wmap(num_nodes);
    for (int u = 0; u < num_nodes; ++u)
        wmap[u].assign(adj[u].size(), 1);

    if (partition.empty()) {
        if (generate_type == 1) {
            std::mt19937 rng(std::random_device{}());
            std::uniform_real_distribution<float> prob(0.0f, 1.0f);
            std::uniform_int_distribution<int> weight(2, 10);
            for (int u = 0; u < num_nodes; ++u)
                for (auto& w : wmap[u])
                    w = (prob(rng) < DEFAULT_RANDOM_WEIGHT_PROBABILITY) ? -weight(rng) : 1;
        }
        // generate_type == 0 → all 1s (already set)
        return wmap;
    }

    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<float> prob(0.0f, 1.0f);

    for (int u = 0; u < num_nodes; ++u) {
        for (int i = 0; i < (int)adj[u].size(); ++i) {
            int v = adj[u][i].to;
            bool is_cut = (partition[u] != partition[v]);
            if (!is_cut) continue;

            int orig = wmap[u][i];
            output += "Edge (" + std::to_string(u) + "," + std::to_string(v) + ") is a cut edge.\n";

            switch (generate_type) {
                case 0:
                    wmap[u][i] = -NEGATIVE_WEIGHT;
                    break;
                case 1:
                    wmap[u][i] = orig - NEGATIVE_WEIGHT;
                    break;
                case 2:
                    wmap[u][i] = -(int)((double)1 / (double)orig * INVERSE_WEIGHT_MULTIPLIER);
                    break;
                case 3:
                    if (prob(rng) < stochastic)
                        wmap[u][i] = -NEGATIVE_WEIGHT;
                    break;
            }
        }
    }
    return wmap;
}

// Top-level function for plain graphs, mirrors generate_weights_from_hgp.
void generate_weights_from_graph(const fs::path& graph_path,
                                  const fs::path& weighted_path,
                                  const fs::path& partition_path,
                                  int generate_type,
                                  bool suppress_output,
                                  double stochastic)
{
    std::string output;
    int num_nodes = 0, num_edges = 0;
    bool had_edge_weights = false;

    AdjacencyList adj = read_metis_graph(graph_path, num_nodes, num_edges, had_edge_weights);

    if (had_edge_weights)
        output += "Warning: input graph already had edge weights; they will be replaced.\n";

    std::vector<PartitionID> partition;
    if (!partition_path.empty())
        partition = readPartitionFile(partition_path, (HypernodeID)num_nodes);

    auto wmap = assign_graph_edge_weights(adj, partition, generate_type, stochastic, output);

    write_metis_graph(weighted_path, adj, wmap);

    if (!suppress_output)
        std::cout << output;
}

// ============================================================
// main
// ============================================================

int main(int argc, char* argv[]) {
    std::string hypergraph_path;
    std::string weighted_dir;
    std::string partition_path;
    double stochastic = 0.5;
    int generate_type = 0;

    CLI::App app;

    app.add_option("-H,--hypergraph", hypergraph_path,
                   "Hypergraph/graph filename or directory")->required();
    app.add_option("-w,--weighted", weighted_dir,
                   "Output directory for weighted files")->required();
    app.add_option("-g,--generate-type", generate_type,
                   "Weight generation mode:\n"
                   "  No partition: 0=all-ones, 1=random\n"
                   "  With partition: 0=cut→-NEG, 1=cut→w-NEG, 2=cut→-1/w*MUL, 3=cut→-NEG(stochastic)");

    CLI::Option* opt       = app.add_option("-p,--partitioned", partition_path,
                                            "Partition file (optional)");
    CLI::Option* opt_stoch = app.add_option("-s,--stochastic", stochastic,
                                            "Probability for stochastic mode (default 0.5)");
    app.add_option("-d,--default-neg", NEGATIVE_WEIGHT,
                   "Penalty for cut edges (default 10)")
       ->default_val(NEGATIVE_WEIGHT)->capture_default_str();
    bool output_print = false;
    app.add_flag("-n,--no-output", output_print,
                 "Suppress per-edge console output")
       ->default_val("false")->capture_default_str();

    CLI11_PARSE(app, argc, argv);

    fs::path weighted_dir_path   = weighted_dir;
    fs::path hypergraph_path_p   = hypergraph_path;
    fs::path partition_path_p    = partition_path;

    if (!fs::is_directory(weighted_dir_path))
        throw fs::filesystem_error("weighted path is not a directory",
                                   weighted_dir_path,
                                   std::make_error_code(std::errc::not_a_directory));

    // Helper lambdas to dispatch by extension
    auto process_hgr = [&](const fs::path& hg, const fs::path& part) {
        fs::path out = weighted_dir_path / hg.filename();
        generate_weights_from_hgp(hg, out, part, generate_type, output_print, stochastic);
    };

    auto process_graph = [&](const fs::path& g, const fs::path& part) {
        fs::path out = weighted_dir_path / g.filename();
        generate_weights_from_graph(g, out, part, generate_type, output_print, stochastic);
    };

    const fs::path part = (opt->count() > 0) ? partition_path_p : fs::path();

    if (fs::is_regular_file(hypergraph_path_p)) {
        auto ext = hypergraph_path_p.extension();
        if (ext == ".hgr") {
            process_hgr(hypergraph_path_p, part);
        } else if (ext == ".graph" || ext == ".metis") {
            process_graph(hypergraph_path_p, part);
        } else {
            std::cerr << "Unsupported file extension: " << ext
                      << "  (expected .hgr, .graph, or .metis)\n";
            return 1;
        }
    } else if (fs::is_directory(hypergraph_path_p)) {
        for (const auto& entry : fs::directory_iterator(hypergraph_path_p)) {
            if (!entry.is_regular_file()) continue;
            auto ext = entry.path().extension();
            if (ext == ".hgr")
                process_hgr(entry.path(), part);
            else if (ext == ".graph" || ext == ".metis")
                process_graph(entry.path(), part);
        }
    } else {
        std::cerr << "Invalid hypergraph path: " << hypergraph_path_p << "\n";
        return 1;
    }

    return 0;
}