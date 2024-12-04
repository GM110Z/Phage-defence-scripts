import pandas as pd
import argparse

# Set up command line argument parsing
parser = argparse.ArgumentParser(description='Process data from PFAM, PADLOC, Antidefense, and AMRFinder.')
parser.add_argument('--pfam', type=str, required=True, help='Path to the PFAM data file (e.g., filtered_pfam_domains_with_desc.tsv)')
parser.add_argument('--padloc', type=str, required=True, help='Path to the PADLOC data file (e.g., allICEspadloc.csv)')
parser.add_argument('--operon', type=str, required=True, help='Path to the operon results file (e.g., operon_results_with_clusters.csv)')
parser.add_argument('--antidefense', type=str, help='Path to the Antidefense data file (e.g., alldefense.txt)')
parser.add_argument('--amrfinder', type=str, help='Path to the AMRFinder data file (e.g., amrall.tsv)')
parser.add_argument('--filter_representative', action='store_true', help='Filter for representative operons')
parser.add_argument('--output', type=str, default='combined_output.txt', help='Output file name (default: combined_output.txt)')

args = parser.parse_args()

# Load the operon results table with clusters
df2 = pd.read_csv(args.operon, sep=",")  # Operon table with Protein and operon info

# Load PFAM domain table
df1 = pd.read_csv(args.pfam, sep="\t")  # PFAM table with domain info

# Load the PADLOC results
df_new = pd.read_csv(args.padloc, sep=",")  # PADLOC output with system and target name

# Initialize merged DataFrame
merged_final_df = pd.DataFrame()

# Merge operon table with the PFAM domains table on Protein_ID and target name
merged_df = pd.merge(df2, df1, left_on='Protein_ID', right_on='target name', how='left')
final_df = merged_df[['Protein_ID', 'nuccore_id', 'start', 'stop', 'strand',
                      'operon_number', 'product', 'accession',
                      'query name', 'E-value', 'description of target']]
merged_final_df = final_df

# Merge with PADLOC data
merged_final_df = pd.merge(merged_final_df, df_new[['system', 'target.name']],
                           left_on='Protein_ID', right_on='target.name',
                           how='left')

# Drop the redundant 'target.name' column
merged_final_df = merged_final_df.drop(columns=['target.name'])

# Merge with Antidefense data if the file is provided
if args.antidefense:
    antidefense_df = pd.read_csv(args.antidefense, sep="\t")  # Antidefense data with sys_beg
    antidefense_df['Antidefense'] = antidefense_df['Antidefense'] + '_' + antidefense_df['AntidefenceSubtype']
    antidefense_df = antidefense_df.drop(columns=['AntidefenceSubtype'])

    merged_final_df = pd.merge(merged_final_df, antidefense_df[['Antidefense', 'sys_beg']],
                               left_on='Protein_ID', right_on='sys_beg',
                               how='left')

    # Drop the redundant 'sys_beg' column
    merged_final_df = merged_final_df.drop(columns=['sys_beg'])

# Merge with AMRFinder data if the file is provided
if args.amrfinder:
    amrfinder_df = pd.read_csv(args.amrfinder, sep='\t')  # AMRFinder data with additional columns
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

# Save the merged data to the specified output file
merged_final_df.to_csv(args.output, index=False, sep='\t')

# Filter for representative operons if the option is selected
if args.filter_representative:
    # Load data from the merged output file
    data = pd.read_csv(args.output, sep='\t')

    # Group by 'MmseqCluster' and filter for representatives
    grouped = data.groupby('MmseqCluster')
    cleaned_data = []

    for mcluster, group in grouped:
        if not group.empty:
            representative = group.iloc[0]  # Select the first entry as representative
            cleaned_data.append(representative)

    # Convert cleaned data back to a DataFrame
    cleaned_df = pd.DataFrame(cleaned_data)

    # Save the cleaned DataFrame
    cleaned_df.to_csv('cleaned_output.csv', sep='\t', index=False)
else:
    print("Representative operon filtering step was not run.")


