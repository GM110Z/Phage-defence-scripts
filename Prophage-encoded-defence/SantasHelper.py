import pandas as pd
import argparse
# Set up command line argument parsing
parser = argparse.ArgumentParser(description='Process data from PFAM, PADLOC, Antidefense, and AMRFinder.')
parser.add_argument('--padloc', action='store_true', help='Include PADLOC data')
parser.add_argument('--pfam', action='store_true', help='Include PFAM data')
parser.add_argument('--antidefense', action='store_true', help='Include Antidefense data')
parser.add_argument('--amrfinder', action='store_true', help='Include AMRFinder data')

args = parser.parse_args()

# Check if at least PFAM and PADLOC are included
if not (args.pfam and args.padloc):
    print("You must include both PADLOC and PFAM data.")
    exit(1)

# Load the PFAM domain table with a description column
df1 = pd.read_csv('filtered_pfam_domains_with_desc.tsv', sep="\t")  # PFAM table with domain info

# Load the operon results table with clusters
df2 = pd.read_csv('operon_results_with_clusters.csv', sep=",")  # Operon table with Protein and operon info

# Load the PADLOC results and reduce to relevant columns
df_new = pd.read_csv('allICEspadloc.csv', sep=",")  # PADLOC output with system and target name

# Initialize merged DataFrame with PFAM and PADLOC data
merged_final_df = pd.DataFrame()

# Merge the operon table with the PFAM domains table on Protein_ID and target name
if args.pfam:
    merged_df = pd.merge(df2, df1, left_on='Protein_ID', right_on='target name', how='left')
    final_df = merged_df[['Protein_ID', 'nuccore_id', 'start', 'stop', 'strand', 
                          'operon_number', 'product', 'accession', 
                          'query name', 'E-value', 'description of target']]
    merged_final_df = final_df

# Merge the resulting DataFrame with the PADLOC data on Protein_ID
if args.padloc:
    merged_final_df = pd.merge(merged_final_df, df_new[['system', 'target.name']], 
                               left_on='Protein_ID', right_on='target.name', 
                               how='left')

    # Drop the redundant 'target.name' column
    merged_final_df = merged_final_df.drop(columns=['target.name'])

# Merge with the updated Antidefense data if the option is selected
if args.antidefense:
    antidefense_df = pd.read_csv('alldefense.txt', sep="\t")  # Antidefense data with sys_beg
    antidefense_df['Antidefense'] = antidefense_df['Antidefense'] + '_' + antidefense_df['AntidefenceSubtype']
    antidefense_df = antidefense_df.drop(columns=['AntidefenceSubtype'])

    merged_final_df = pd.merge(merged_final_df, antidefense_df[['Antidefense', 'sys_beg']],
                               left_on='Protein_ID', right_on='sys_beg', 
                               how='left')

    # Drop the redundant 'sys_beg' column from the final DataFrame
    merged_final_df = merged_final_df.drop(columns=['sys_beg'])

# Merge with AMRFinder data if the option is selected
if args.amrfinder:
    amrfinder_df = pd.read_csv('amrall.tsv', sep='\t')  # AMRFinder data with additional columns
    amrfinder_df['AMRFinder'] = amrfinder_df['Sequence name'] + '_' + \
                                amrfinder_df['Scope'] + '_' + \
                                amrfinder_df['Element type'] + '_' + \
                                amrfinder_df['Element subtype'] + '_' + \
                                amrfinder_df['Class'] + '_' + \
                                amrfinder_df['Subclass']
    amrfinder_df = amrfinder_df[['Protein identifier', 'AMRFinder']]

    merged_final_df = pd.merge(merged_final_df, amrfinder_df, 
                               left_on='Protein_ID', right_on='Protein identifier', 
                               how='left')

    # Drop the redundant 'Protein identifier' column
    merged_final_df = merged_final_df.drop(columns=['Protein identifier'])

# Save the final merged data to a new CSV file
merged_final_df.to_csv('combined_with_amrfinder.txt', index=False)

# Load your data
data = pd.read_csv('combined_with_amrfinder.txt', sep='\t')

# Optionally filter for representative operons
if args.filter_representative:
    # Step 1: Group by 'MmseqCluster'
    grouped = data.groupby('MmseqCluster')
    cleaned_data = []

    # Step 2: Iterate over each group
    for mcluster, group in grouped:
        # Step 3: Select the subset for the current cluster
        mcluster_subset = group[group['MmseqCluster'] == mcluster]

        # Step 4: Check if the subset is empty
        if not mcluster_subset.empty:
            # Keep only the first entry or based on your preferred logic
            representative = mcluster_subset.iloc[0]
            cleaned_data.append(representative)

    # Convert cleaned data back to DataFrame
    cleaned_df = pd.DataFrame(cleaned_data)

    # Save or display the cleaned DataFrame
    cleaned_df.to_csv('cleaned_output.csv', sep='\t', index=False)
else:
    print("Representative operon filtering step was not run.")



