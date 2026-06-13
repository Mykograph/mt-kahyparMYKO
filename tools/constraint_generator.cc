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

using namespace mt_kahypar;
namespace fs = std::filesystem;

const std::string CONSTRAINT_FILE_EXTENSION = ".constraints.txt";

inline std::string to_trimmed_string(double x) {
    std::ostringstream oss;
    oss << x;
    std::string s = oss.str();

    if (s.find('.') != std::string::npos) {
        s.erase(s.find_last_not_of('0') + 1);
        if (s.back() == '.') {
            s.pop_back();
        }
    }
    return s;
}

std::optional<fs::path> findFileWithPrefix(const fs::path& path, const std::string& prefix) {
    if (fs::is_regular_file(path)) {
        auto name = path.filename().string();
        if (name.rfind(prefix, 0) == 0) {  // starts with prefix
            return path;
        }
    }
    if (fs::is_directory(path)) {
        for (const auto& entry : fs::directory_iterator(path)) {
            if (entry.is_regular_file()) {
                auto filename = entry.path().filename().string();
                if (filename.rfind(prefix, 0) == 0) {
                    return entry.path();
                }
            }
        }
    }
    return std::nullopt;
}

void print_constraints(fs::path constraints_path, 
                        vec<std::pair<HypernodeID, HypernodeID>>& constraint_list, 
                        HypernodeID num_constraints) {
    std::mt19937 gen(12345);
    std::shuffle(constraint_list.begin(), constraint_list.end(), gen);

    std::ofstream out_stream(constraints_path);
    for (HypernodeID i = 0; i < num_constraints; i++) {
        std::pair<HypernodeID, HypernodeID> constraint = constraint_list[i];
        out_stream << constraint.first << " " << constraint.second << std::endl;
    }
    out_stream.close();
}

HypernodeID pick_random(const HypernodeID limit) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist(0, limit);
    return dist(gen);
}

bool try_add_constraint(HypernodeID a,
                        HypernodeID b,
                        const HypernodeID max_constraints_per_node,
                        std::unordered_map<HypernodeID, HypernodeID>& constraints_count_per_node,
                        std::unordered_map<HypernodeID, std::unordered_set<HypernodeID>>& constraints_per_node) {
    if (a == b) return false;
    if (constraints_count_per_node[a] >= max_constraints_per_node) return false;
    if (constraints_count_per_node[b] >= max_constraints_per_node) return false;

    auto& set_a = constraints_per_node[a];
    auto& set_b = constraints_per_node[b];
    if (set_a.find(b) != set_a.end() || set_b.find(a) != set_b.end()) return false;

    set_a.insert(b);
    set_b.insert(a);
    constraints_count_per_node[a]++;
    constraints_count_per_node[b]++;
    return true;
}

HypernodeID generate_constraints_from_hg(const fs::path hg_path,
                                        const fs::path constraints_path,
                                        const float constraints_percentage,
                                        const HypernodeID max_constraints_per_node,
                                        const HypernodeID desired_node_degree) {
    // Read Hypergraph
    HyperedgeID num_edges;
    HypernodeID num_nodes;
    HyperedgeID num_removed_single_pin_hyperedges;
    io::HyperedgeVector hyperedges;
    vec<HyperedgeWeight> hyperedges_weight;
    vec<HypernodeWeight> hypernodes_weight;

    io::readHypergraphFile(hg_path.string(), num_edges, num_nodes, num_removed_single_pin_hyperedges,
                         hyperedges, hyperedges_weight, hypernodes_weight, false, false);
    ALWAYS_ASSERT(hyperedges.size() == num_edges);
    ALWAYS_ASSERT(num_removed_single_pin_hyperedges == 0);

    HypernodeID num_constraints = num_nodes * constraints_percentage;
    HypernodeID constraint_count = 0;
    std::unordered_map<HypernodeID, HypernodeID> constraints_count_per_node;
    std::unordered_map<HypernodeID, std::unordered_set<HypernodeID>> constraints_per_node;

    std::ofstream out_stream(constraints_path);
    while(constraint_count < num_constraints) {
        HypernodeID node = pick_random(num_nodes - 1);
        HypernodeID node_constraints = ((desired_node_degree == 0)? pick_random(max_constraints_per_node) : desired_node_degree);
        HypernodeID count_node_constraints = 0;
        while (count_node_constraints < node_constraints && constraint_count < num_constraints) {
            HypernodeID other_node = pick_random(num_nodes - 1);
            if (constraints_count_per_node[node] >= max_constraints_per_node) break;
            if (try_add_constraint(
                node,
                other_node,
                max_constraints_per_node,
                constraints_count_per_node,
                constraints_per_node
            )) {
                count_node_constraints++;
                constraint_count++;
                out_stream << node << " " << other_node << std::endl;
            }
        }
    }
    out_stream.close();
    return constraint_count;
}

HypernodeID generate_constraints_from_partitioned_hg(const fs::path hg_path, 
                                                        const fs::path part_hg_path, 
                                                        const fs::path constraints_path, 
                                                        const float constraints_percentage, 
                                                        const HypernodeID max_constraints_per_node,
                                                        const HypernodeID desired_node_degree) {
    // Read Hypergraph header and partitioned Hypergraph
    HyperedgeID num_edges = 0;
    HypernodeID num_nodes = 0;
    HyperedgeID num_removed_single_pin_hyperedges_1 = 0;
    io::HyperedgeVector hyperedges_1;
    vec<HyperedgeWeight> hyperedges_weight_1;
    vec<HypernodeWeight> hypernodes_weight_1;
    std::vector<PartitionID> partitions;
    io::readHypergraphFile(hg_path.string(), num_edges, num_nodes, num_removed_single_pin_hyperedges_1,
                         hyperedges_1, hyperedges_weight_1, hypernodes_weight_1, false, false);
    io::readPartitionFile(part_hg_path.string(), num_nodes, partitions);

    HypernodeID num_constraints = num_nodes * constraints_percentage;
    HypernodeID constraint_count = 0;
    std::unordered_map<HypernodeID, HypernodeID> constraints_count_per_node;
    std::unordered_map<HypernodeID, std::unordered_set<HypernodeID>> constraints_per_node;

    std::ofstream out_stream(constraints_path);
    while(constraint_count < num_constraints) {
        HypernodeID node = pick_random(num_nodes - 1);
        HypernodeID node_constraints = ((desired_node_degree == 0)? pick_random(max_constraints_per_node) : desired_node_degree);
        PartitionID node_partition = partitions[node];
        HypernodeID count_node_constraints = 0;
        while (count_node_constraints < node_constraints && constraint_count < num_constraints) {
            HypernodeID other_node = pick_random(num_nodes - 1);
            if (constraints_count_per_node[node] >= max_constraints_per_node) break;
            if (node_partition != partitions[other_node] &&
                try_add_constraint(
                node,
                other_node,
                max_constraints_per_node,
                constraints_count_per_node,
                constraints_per_node
            )) {
                count_node_constraints++;
                constraint_count++;
                out_stream << node << " " << other_node << std::endl;
            }
        }
    }
    out_stream.close();
    return constraint_count;
}

HypernodeID generate_constraints_from_partitioned_hg_old_way(const fs::path hg_path, 
                                                        const fs::path part_hg_path, 
                                                        const fs::path constraints_path, 
                                                        const float constraints_percentage, 
                                                        const HypernodeID max_constraints_per_node,
                                                        const HypernodeID desired_node_degree) {
    unused(desired_node_degree);
    HyperedgeID num_edges = 0;
    HypernodeID num_nodes = 0;
    HyperedgeID num_removed_single_pin_hyperedges_2 = 0;
    io::HyperedgeVector hyperedges_2;
    vec<HyperedgeWeight> hyperedges_weight_2;
    vec<HypernodeWeight> hypernodes_weight_2;
    std::vector<PartitionID> partitions;
    io::readHypergraphFile(hg_path.string(), num_edges, num_nodes, num_removed_single_pin_hyperedges_2,
                         hyperedges_2, hyperedges_weight_2, hypernodes_weight_2, false, false);
    io::readPartitionFile(part_hg_path.string(), num_nodes, partitions);

    HypernodeID num_constraints = num_nodes * constraints_percentage;

    vec<std::pair<HypernodeID, HypernodeID>> constraint_list;
    constraint_list.reserve(num_nodes * max_constraints_per_node);
    HypernodeID constraint_count = 0;
    std::unordered_map<HypernodeID, HypernodeID> constraints_per_node;
    for (HypernodeID node = 0; node < num_nodes; node++) {
        PartitionID node_partition = partitions[node];
        for (HypernodeID other_node = node + 1; other_node < num_nodes; other_node++) {
            if (constraints_per_node[node] >= max_constraints_per_node) break;
            if (constraints_per_node[other_node] >= max_constraints_per_node) continue;

            if (node_partition != partitions[other_node]) {
                constraint_count++;
                constraints_per_node[node]++;
                constraints_per_node[other_node]++;
                constraint_list.push_back(std::make_pair(node, other_node));
            }
        }
        constraints_per_node.erase(node);
    }
    num_constraints = std::min(num_constraints, constraint_count);
    print_constraints(constraints_path, constraint_list, num_constraints);
    return num_constraints;
}

int main(int argc, char* argv[]) {
    fs::path hypergraph_path;
    vec<fs::path> hypergraph_files;
    std::optional<fs::path> partitioned_hypergraph_path;
    fs::path constraint_dir;
    HypernodeID k;
    HypernodeID desired_node_degree;
    float num_constraints_percentage;
    bool generate_from_partitioned_hg_old_way;
    HypernodeID max_constraints_per_node;

    std::string hypergraph_path_str;
    std::string constraint_dir_str;
    std::string part_hg_path_str;

    CLI::App app{"Constraint Generator"};
    app.add_option("-g,--hypergraph", hypergraph_path_str, "Hypergraph Filename or directory")->required()->type_name("<path>");
    app.add_option("-c,--constraint", constraint_dir_str, "Constraint directory")->required()->type_name("<path>");
    app.add_option("-k,--blocks", k, "Number of blocks")->required()->type_name("<int>");
    app.add_option("-d,--degree", desired_node_degree, "Desired constraint node degree")->default_val(0)->type_name("<int>");
    app.add_option("-n,--num-constraints", num_constraints_percentage, "Number of constraints (optional)")->default_val(1.0)->type_name("<float>");
    app.add_flag("--old", generate_from_partitioned_hg_old_way, "Generate constraints from part hg the old way (optional)");
    app.add_option("-p,--part-hypergraph", part_hg_path_str, "Partitioned-Hypergraph Filename or directory (optional)")->type_name("<path>");

    CLI11_PARSE(app, argc, argv);

    hypergraph_path = hypergraph_path_str;
    constraint_dir = constraint_dir_str;
    if (!part_hg_path_str.empty()) {
        partitioned_hypergraph_path = fs::path(part_hg_path_str);
    }

    max_constraints_per_node = k - 1;

    if (!fs::is_directory(constraint_dir)) {
        throw fs::filesystem_error(
            "constraint path is not a directory",
            constraint_dir,
            std::make_error_code(std::errc::not_a_directory)
        );
    }

    if (fs::is_regular_file(hypergraph_path) && hypergraph_path.extension() == ".hgr") {
        hypergraph_files.push_back(hypergraph_path);
    } else if (fs::is_directory(hypergraph_path)) {
        for (const auto& file : fs::directory_iterator(hypergraph_path)) {
            if (fs::is_regular_file(file.path()) && file.path().extension() == ".hgr") {
                hypergraph_files.push_back(file.path());
            }
        }
    } else {
        throw fs::filesystem_error(
            "hypergraph file doesnt exist",
            hypergraph_path,
            std::make_error_code(std::errc::not_a_directory)
        );
    }

    for (const auto& hg_file : hypergraph_files) {
        fs::path constraint_file = constraint_dir / (hg_file.filename().string() + "." + std::to_string(k) + ".txt");
        HypernodeID generated_constraints;
        if (partitioned_hypergraph_path) {
            auto part_hg_file = findFileWithPrefix(partitioned_hypergraph_path.value(), hg_file.filename().string() + ".part" + std::to_string(k));
            if (!part_hg_file){
                throw fs::filesystem_error(
                    "no matching partitioned hypergraph file found in ",
                    partitioned_hypergraph_path.value(),
                    std::make_error_code(std::errc::not_a_directory)
                );
            }
            if(generate_from_partitioned_hg_old_way) {
                constraint_file.replace_filename(constraint_file.stem().string() + ".constraints" + constraint_file.extension().string());
                generated_constraints = generate_constraints_from_partitioned_hg_old_way(hg_file, part_hg_file.value(), constraint_file, num_constraints_percentage, max_constraints_per_node, desired_node_degree);
            } else {
                constraint_file.replace_filename(constraint_file.stem().string() + ".constraints" + constraint_file.extension().string());
                generated_constraints = generate_constraints_from_partitioned_hg(hg_file, part_hg_file.value(), constraint_file, num_constraints_percentage, max_constraints_per_node, desired_node_degree);
            }
        } else {
            constraint_file.replace_filename(constraint_file.stem().string() + ".constraints" + constraint_file.extension().string());
            generated_constraints = generate_constraints_from_hg(hg_file, constraint_file, num_constraints_percentage, max_constraints_per_node, desired_node_degree);
        }
        LOG << "";
        LOG << "Generated " << generated_constraints << "constraints for hg:" << hg_file.filename();
    }
    LOG << "";
    return 0;
}