import pandas as pd

table1=pd.read_csv("All-pros.tsv", sep="\t")
table2=pd.read_csv("defence-prophage-clusters.tsv", sep="\t")
# Rename columns to make sure they are correctly aligned for operations, if needed
table1.columns = ['genome', 'cluster_id']
table2.columns = ['genome']

# Merging tables on the 'genome' to find matches
initial_matches = pd.merge(table1, table2, on='genome')

# Get unique cluster IDs from matches to select all related entries from table1
clusters_to_select = initial_matches['cluster_id'].unique()

# Select all entries from table1 that have a cluster_id matching any of the matched cluster_ids
final_result = table1[table1['cluster_id'].isin(clusters_to_select)]

# Saving the result to a new CSV file
final_result.to_csv('expanded_matches.csv', index=False)

