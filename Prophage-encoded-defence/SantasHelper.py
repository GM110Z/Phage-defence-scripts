#combine operons predictions with pfams and padloc
import pandas as pd

# Assuming the files are tab-separated for both tables
df1 = pd.read_csv('filtered_pfam_domains_with_desc.tsv', sep="\t")  # The PFAM table with Domain information
df2 = pd.read_csv('representative_operons.csv', sep=",")  # The second table with Protein and operon information
df_new = pd.read_csv('padloc-short.txt', sep="\t")  # Combine PADLOC output of all files and reduce columns to just 'system' and 'target.name'

# Perform the merge based on 'Protein_ID' from df2 and 'Domain_ID' from df1
merged_df = pd.merge(df2, df1,  left_on='Protein_ID', right_on='Domain_ID', how='left')

# Now select only the columns from df2 and specific columns from df1 that you want
final_df = merged_df[['Protein_ID', 'nuccore_id', 'start', 'stop', 'strand', 'operon_number', 'Cluster_ID', 'Description_of_Target', 'E_Value']]
#now add padloc data

# Merge with the new system table based on matching 'Protein_ID' and 'target.name'
final_df = pd.merge(final_df, df_new[['system', 'target.name']], left_on='Protein_ID', right_on='target.name', how='left')# Drop duplicates based on all columns (if values in all these columns are the same) to eliminate duplications of the same nuccore
final_df = final_df.drop_duplicates(subset=['Protein_ID', 'nuccore_id', 'start', 'stop', 'strand', 'operon_number', 'Cluster_ID', 'Description_of_Target', 'E_Value'])
final_df = final_df.drop_duplicates(subset=['Protein_ID', 'nuccore_id', 'start', 'stop', 'strand', 'operon_number', 'Cluster_ID', 'Description_of_Target'])

#drop duplicates of the same operon on different nuccore
final_df = final_df.drop_duplicates(subset=['Protein_ID', 'Cluster_ID', 'Description_of_Target', 'E_Value'])


# Save the result to a new file or display it
final_df.to_csv('Knonw-and-New-candidates-in-hotspot.csv', index=False)

