#!/usr/bin/env python3
import os
import csv

BASE_DIR    = os.path.expanduser("~/mt-kahypar/generated_instances")
RESULTS_DIR = os.path.expanduser("~/mt-kahypar/generated_instances")
k_values    = [2, 16]
percentages = range(10, 101, 10)
neg_values  = [1, 10, 100]
ALGORITHM   = "Mt-KaHyPar-D"
EPSILON     = 0.03
THREADS     = 8

def parse_hmetis(filepath):
    hyperedges   = []
    edge_weights = []
    with open(filepath) as f:
        lines = [l for l in f if not l.startswith('%')]
    header           = lines[0].split()
    num_edges        = int(header[0])
    num_vertices     = int(header[1])
    has_edge_weights = len(header) > 2 and int(header[2]) in (1, 11)
    for line in lines[1:num_edges+1]:
        parts = list(map(int, line.split()))
        if has_edge_weights:
            edge_weights.append(parts[0])
            parts = parts[1:]
        else:
            edge_weights.append(1)
        hyperedges.append([v - 1 for v in parts])
    return hyperedges, edge_weights, num_vertices

def parse_partition(filepath):
    with open(filepath) as f:
        return [int(line.strip()) for line in f if line.strip()]

def compute_km1(hyperedges, edge_weights, partition):
    km1 = 0
    for edge, w in zip(hyperedges, edge_weights):
        blocks = set(partition[v] for v in edge)
        km1 += (len(blocks) - 1) * w
    return km1

def compute_cut(hyperedges, edge_weights, partition):
    cut = 0
    for edge, w in zip(hyperedges, edge_weights):
        blocks = set(partition[v] for v in edge)
        if len(blocks) > 1:
            cut += w
    return cut

def count_negative_cut_edges(hyperedges, edge_weights, partition):
    neg_cut_count  = 0
    neg_cut_weight = 0
    for edge, w in zip(hyperedges, edge_weights):
        if w >= 0:
            continue
        blocks = set(partition[v] for v in edge)
        if len(blocks) > 1:
            neg_cut_count  += 1
            neg_cut_weight += w
    return neg_cut_count, neg_cut_weight

def count_positive_cut_edges(hyperedges, edge_weights, partition):
    pos_cut_count  = 0
    pos_cut_weight = 0
    for edge, w in zip(hyperedges, edge_weights):
        if w <= 0:
            continue
        blocks = set(partition[v] for v in edge)
        if len(blocks) > 1:
            pos_cut_count  += 1
            pos_cut_weight += w
    return pos_cut_count, pos_cut_weight

def find_partition_file(partition_dir, graph_name):
    direct = os.path.join(partition_dir, graph_name)
    if os.path.exists(direct):
        return direct
    base = graph_name[:graph_name.index('.seed')] if '.seed' in graph_name else graph_name
    try:
        for fname in os.listdir(partition_dir):
            if fname.startswith(base):
                return os.path.join(partition_dir, fname)
    except FileNotFoundError:
        pass
    return None

def process(hgr_dir, partition_dir, output_file, k):
    if not os.path.exists(hgr_dir):
        print(f"  Skipping (not found): {hgr_dir}")
        return
    if not os.path.exists(partition_dir):
        print(f"  Skipping (not found): {partition_dir}")
        return
    rows = []
    for graph_name in sorted(os.listdir(hgr_dir)):
        hgr_file = os.path.join(hgr_dir, graph_name)
        if not os.path.isfile(hgr_file):
            continue
        part_file = find_partition_file(partition_dir, graph_name)
        if part_file is None:
            print(f"  Warning: Partition not found for {graph_name} in {partition_dir}")
            continue
        try:
            hyperedges, edge_weights, num_vertices = parse_hmetis(hgr_file)
            partition = parse_partition(part_file)
            if len(partition) != num_vertices:
                print(f"  Warning: vertex count mismatch for {graph_name}")
                continue
            km1                           = compute_km1(hyperedges, edge_weights, partition)
            cut                           = compute_cut(hyperedges, edge_weights, partition)
            neg_cut_count, neg_cut_weight = count_negative_cut_edges(hyperedges, edge_weights, partition)
            pos_cut_count, pos_cut_weight = count_positive_cut_edges(hyperedges, edge_weights, partition)
            rows.append([
                ALGORITHM, graph_name, "no", 0, k, EPSILON, THREADS,
                0.0, 0, "km1", km1, cut, "no",
                neg_cut_count, neg_cut_weight,
                pos_cut_count, pos_cut_weight
            ])
        except Exception as e:
            print(f"  Error processing {graph_name}: {e}")
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(rows)
    print(f"  Written: {output_file}")

for k in k_values:
    for p in percentages:
        for neg in neg_values:
            hgr_dir    = os.path.join(BASE_DIR, f"hg-{k}-{p}-{neg}")
            part_dir_a = os.path.join(BASE_DIR, f"partition-{k}")
            part_dir_b = os.path.join(BASE_DIR, f"partition-{k}-B")
            output_a   = os.path.join(RESULTS_DIR, f"initial-{k}-{p}-{neg}.csv")
            output_b   = os.path.join(RESULTS_DIR, f"baseline-{k}-{p}-{neg}.csv")
            print(f"Processing k={k} p={p} neg={neg}")
            process(hgr_dir, part_dir_a, output_a, k)
            process(hgr_dir, part_dir_b, output_b, k)
