#!/usr/bin/env python


import sys
import pandas as pd

# Get argument
args = sys.argv
matrix_tsv_file = args[1]  # ANIclustermap_matrix.tsv
cluster_ani_thr = float(args[2])  # e.g. 95.0
cluster_tsv_file = args[3]  # Output cluster table file

# Parse cluster ANI matrix
df = pd.read_table(matrix_tsv_file)
cluster_id = 1
cluster_base_idx = 0
cluster_size_record = 1
genome_name2cluster_id = {}
for i, genome_name in enumerate(df.columns):
    cluster_candidate_df = df.iloc[cluster_base_idx : i + 1, cluster_base_idx : i + 1]
    ani_thr_match_count = (cluster_candidate_df > cluster_ani_thr).sum().sum()
    if ani_thr_match_count != cluster_size_record**2:
        cluster_id += 1
        cluster_base_idx = i
        cluster_size_record = 1

    genome_name2cluster_id[genome_name] = cluster_id
    cluster_size_record += 1

# Output cluster table
cluster_table_dict = {
    "genome": genome_name2cluster_id.keys(),
    "cluster_id": genome_name2cluster_id.values(),
}
cluster_table_df = pd.DataFrame(cluster_table_dict)
cluster_table_df.to_csv(cluster_tsv_file, sep="\t", index=False)
