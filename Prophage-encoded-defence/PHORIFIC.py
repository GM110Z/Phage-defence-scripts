#phorific-final

import sys
from Bio import SeqIO
import csv
import glob
import os
from collections import defaultdict
import pandas as pd
import numpy as np

## Step 1: Edit MMseqs clustering table to introduce cluster numbers
def load_clustering_table(file_path):
    return pd.read_csv(file_path, sep="\t", header=None, names=["ClusterRep", "Protein_ID"])

def assign_numerical_ids(df):
    unique_clusters = df['ClusterRep'].unique()
    cluster_to_id = {cluster: idx + 1 for idx, cluster in enumerate(unique_clusters)}
    df['ClusterRep'] = df['ClusterRep'].map(cluster_to_id)
    return df, cluster_to_id

def save_updated_table(df, file_path):
    df.to_csv(file_path, sep="\t", index=False, header=False)

# Load MMseqs cluster data to map Protein_ID to ClusterRep
def load_cluster_data(cluster_file):
    cluster_df = pd.read_csv(cluster_file, sep="\t", header=None, names=["ClusterRep", "Protein_ID"])
    return cluster_df.set_index("Protein_ID")["ClusterRep"].to_dict()

# Part 2: Genbank files to table of features
def combine_genbank_files():
    for f in glob.glob("*.gb"):
        os.system("cat " + f + " >> OutFile.gb")

def parse_combined_records():
    combined_files = SeqIO.parse("OutFile.gb", format="genbank")
    with open('prefiltered.txt', mode='w') as parsed_output:
        parsed_output = csv.writer(parsed_output, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        parsed_output.writerow(['Protein_ID', 'product', 'molecule', 'start', 'end', 'strand'])
        for rec in combined_files:
            contig = rec.id
            for feature in rec.features:
                if feature.type == 'CDS':
                    try:
                        protein_id = feature.qualifiers['protein_id'][0].strip()  # Strip whitespace
                        product = feature.qualifiers['product'][0].strip() if 'product' in feature.qualifiers else 'Unknown'
                        print(protein_id, "\t", product, "\t", contig, "\t", feature.location.start, "\t", feature.location.end, "\t", feature.location.strand, "\n", file=open('prefiltered.txt', 'a'))
                    except KeyError:
                        print("Pseudogene", "\t", "Unknown", "\t", contig, "\t", feature.location.start, "\t", feature.location.end, "\t", feature.location.strand, "\n", file=open('prefiltered.txt', 'a'))

def clean_prefiltered_file():
    df = pd.read_csv('prefiltered.txt', sep="\t", low_memory=False, skipinitialspace=True)

    # Ensure 'start' and 'end' columns are strings
    df['start'] = df['start'].astype(str)
    df['end'] = df['end'].astype(str)

    # Proceed with filtering based on string contains
    df2 = df[~df['start'].str.contains('>')]
    df3 = df2[~df2['end'].str.contains('>')]
    df4 = df3[~df3['start'].str.contains('<')]
    df5 = df4[~df4['end'].str.contains('<')]

    # Adjust strand and molecule formatting
    df5['strand'] = np.where(df5['strand'] == -1, 0, df5['strand'])
    df5['molecule'] = df5['molecule'].str.replace(' ', '')

    # Save the cleaned file
    df5.to_csv("TableOfFeatures.txt", sep="\t", index=False, header=True)

def load_genomic_coordinates(filename):
    if not os.path.exists(filename):
        print(f"Error: {filename} not found.")
        return []
    
    genomic_data = []
    with open(filename, 'r') as file:
        reader = csv.DictReader(file, delimiter="\t")
        for row in reader:
            # Strip whitespace from Protein_ID
            protein_id = row['Protein_ID'].strip()
            # Ensure that start and end are valid integers
            if row['start'] and row['end'] and pd.notna(row['start']) and pd.notna(row['end']):
                try:
                    start = int(row['start'])
                    end = int(row['end'])
                except ValueError:
                    print(f"Skipping row due to invalid start or end values: {row}")
                    continue

                genomic_data.append({
                    'Protein_ID': protein_id,  # Use stripped Protein_ID
                    'start': start,
                    'stop': end,
                    'strand': row['strand'],
                    'nuccore_id': row['molecule'],
                    "product": row['product']
                })
            else:
                print(f"Skipping row due to missing fields: {row}")
    return genomic_data

# Identify operons based on proximity and strand compatibility
def identify_operons(genomic_data, cluster_data, threshold=35):
    operon_assignments = []
    operon_number = 0

    # Group by nuccore_id
    grouped_by_nuccore = defaultdict(list)
    for protein in genomic_data:
        grouped_by_nuccore[protein['nuccore_id']].append(protein)

    for nuccore_id, proteins in grouped_by_nuccore.items():
        sorted_proteins = sorted(proteins, key=lambda x: x['start'])
        current_operon = []
        for i in range(len(sorted_proteins)):
            if not current_operon:
                current_operon.append(sorted_proteins[i])
            else:
                last_protein = current_operon[-1]
                if (sorted_proteins[i]['strand'] == last_protein['strand'] and
                        (sorted_proteins[i]['start'] - last_protein['stop']) <= threshold):
                    current_operon.append(sorted_proteins[i])
                else:
                    operon_number += 1
                    for protein in current_operon:
                        operon_assignments.append({
                            'Protein_ID': protein['Protein_ID'],
                            'nuccore_id': nuccore_id,
                            'start': protein['start'],
                            'stop': protein['stop'],
                            'strand': protein['strand'],
                            'operon_number': operon_number,
                            'product': protein['product'],
                            'ClusterRep': cluster_data.get(protein['Protein_ID'], "Unmatched")  # Improved mapping check
                        })
                    current_operon = [sorted_proteins[i]]

        if current_operon:
            operon_number += 1
            for protein in current_operon:
                operon_assignments.append({
                    'Protein_ID': protein['Protein_ID'],
                    'nuccore_id': nuccore_id,
                    'start': protein['start'],
                    'stop': protein['stop'],
                    'strand': protein['strand'],
                    'operon_number': operon_number,
                    'product': protein['product'],
                    'ClusterRep': cluster_data.get(protein['Protein_ID'], "Unmatched")  # Improved mapping check
                })

    return operon_assignments

# Write operon results to a file
def write_results(output_file, operon_assignments):
    if not operon_assignments:
        print("No operon data to write.")
        return
    with open(output_file, 'w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=['Protein_ID', 'nuccore_id', 'start', 'stop', 'strand', 'product', 'operon_number', 'ClusterRep'])
        writer.writeheader()
        writer.writerows(operon_assignments)

# Main execution flow
if __name__ == "__main__":
    # Paths to input and output files
    original_file_path = 'MMseqs-cluster.tsv'  # path to your original MMseqs clustering file
    updated_file_path = 'MMseqs-clusternumb.tsv'  # path to save updated table with cluster numbers

    # Step 1: Load and map cluster data
    clustering_df = load_clustering_table(original_file_path)
    updated_df, _ = assign_numerical_ids(clustering_df)
    save_updated_table(updated_df, updated_file_path)
    cluster_data = load_cluster_data(updated_file_path)

    # Ensure genbank files are combined and parsed into TableOfFeatures
    combine_genbank_files()
    parse_combined_records()
    clean_prefiltered_file()

    # Part 2: Load genomic data and identify operons with cluster information
    input_file = 'TableOfFeatures.txt'
    genomic_data = load_genomic_coordinates(input_file)
    if not genomic_data:
        print("Genomic data could not be loaded.")
        sys.exit(1)

    operon_assignments = identify_operons(genomic_data, cluster_data)

    # Step 3: Write full operon results
    output_file = 'operon_results_with_clusters.csv'
    write_results(output_file, operon_assignments)

    print("Operon results with clusters saved.")
# Cleanup temporary files
os.system("rm OutFile.gb") 
os.system("rm prefiltered.txt")


