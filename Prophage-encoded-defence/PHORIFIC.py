import sys
from Bio import SeqIO, SeqFeature
import csv
import glob
import os
from collections import defaultdict
import pandas as pd
import numpy as np

##Step 1: Edit MMseqs clustering table to introduce cluster numbers
def load_clustering_table(file_path):
    return pd.read_csv(file_path, sep="\t", header=None, names=["ClusterRep", "SingleMember"])

def assign_numerical_ids(df):
    unique_clusters = df['ClusterRep'].unique()
    cluster_to_id = {cluster: idx + 1 for idx, cluster in enumerate(unique_clusters)}
    df['ClusterRep'] = df['ClusterRep'].map(cluster_to_id)
    return df, cluster_to_id

def save_updated_table(df, file_path):
    df.to_csv(file_path, sep="\t", index=False, header=False)
# Main workflow
if __name__ == "__main__":
# Paths to files
    original_file_path = '/Users/giusym/Desktop/Pseudomonas-projects/Prophages/workdir/clusterRes_cluster.tsv' #edit to your path to input Mmseqs2 table
    updated_file_path = '/Users/giusym/Desktop/Pseudomonas-projects/Prophages/workdir/MMseqs-clusternumb.tsv' #edit to your path for the file you want to create
# Step 1: Edit MMseqs clustering table
clustering_df = load_clustering_table(original_file_path)
updated_df, _ = assign_numerical_ids(clustering_df)
save_updated_table(updated_df, updated_file_path)

#P#art 2: Genbank files to table of features
print('Working on the genomes...')
#combine gbk file in working dir
for f in glob.glob("*.gb"):
         os.system("cat "+f+" >> OutFile.gb")

#Parse combined records
combined_files = SeqIO.parse("OutFile.gb", format="genbank")

with open('prefiltered.txt', mode='w') as parsed_output:   #write table for Suppl.data with hitsID, scores and lenght
                parsed_output = csv.writer(parsed_output, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)
                parsed_output.writerow(['Protein_ID','product', 'molecule', 'start', 'end', 'strand'])
# rint the info of each record
for rec in combined_files:
    contig=rec.id
    for feature in rec.features:
        if feature.type == 'CDS':
            try:
                print(feature.qualifiers['protein_id'][0],"\t",feature.qualifiers['product'][0],"\t",contig,"\t",feature.location.start,"\t",feature.location.end,"\t",feature.location.strand,"\n", file=open('prefiltered.txt', 'a'))
            except KeyError:
                print("Pseudogene","\t","Unknown","\t",contig,"\t",feature.location.start,"\t",feature.location.end,"\t",feature.location.strand,"\n",file=open('prefiltered.txt', 'a'))

df = pd.read_csv('prefiltered.txt', sep="\t", low_memory=False,skipinitialspace=True)
df2=df[~df['start'].str.contains('>')]
df3=df2[~df2['end'].str.contains('>')]
df4=df3[~df3['start'].str.contains('<')]
df5=df4[~df4['start'].str.contains('<')]
df5['strand'] = np.where(df5['strand']== -1, 0, df5['strand'])
df5['molecule'].str.replace(' ', '')
df5.to_csv("TableOfFeatures.txt", "\t", index=False, header=True)

##Part 3: eliminate the hotspot boundaries 
df = pd.read_csv('TableOfFeatures.txt', sep="\t")

# Function to filter based on protein ID boundaries
def eliminate_boundaries_by_protein_id(df):
    def is_boundary(row):
        # Check if the current row is the start or end position for the same protein ID within the same molecule
        protein_rows = df[df['Protein_ID'] == row['Protein_ID']]

        # Get the start and end positions for this protein ID
        protein_start = protein_rows['start'].min()
        protein_end = protein_rows['end'].max()

        # Check if the current row is at the boundary
        return row['start'] == protein_start or row['end'] == protein_end

    # Apply the check and filter out the rows that match
    return df[~df.apply(is_boundary, axis=1)]

# Apply the function to eliminate the desired rows
filtered_df = eliminate_boundaries_by_protein_id(df)

# Save the filtered DataFrame to a new CSV file
filtered_df.to_csv('filtered_file_no_boundaries.csv', index=False)
 
#removing combined file                
os.system("rm OutFile.gb") 
os.system("rm prefiltered.txt") 

# Load genomic coordinates from a CSV file
def load_genomic_coordinates(filename):
    genomic_data = []
    try:
        with open(filename, 'r') as file:
            reader = csv.DictReader(file)
            for row in reader:
                # Ensure all necessary fields are present
                if 'Protein_ID' in row and 'start' in row and 'end' in row and 'strand' in row and 'molecule' in row:
                    genomic_data.append({
                        'Protein_ID': row['Protein_ID'],
                        'start': int(row['start']),
                        'stop': int(row['end']),  # Adjusting to 'end'
                        'strand': row['strand'],
                        'nuccore_id': row['molecule']  # Adjusting to 'molecule'
                    })
                else:
                    print("Missing fields in row:", row)
    except Exception as e:
        print("Error reading the file:", e)
    return genomic_data

# Determine operons based on proximity and strand compatibility
def identify_operons(genomic_data, threshold=35): #change the threshold to the maximum distance you want allow between genes of the same operon.
    operon_assignments = []
    operon_number = 0

    # Group by nuccore_id
    grouped_by_nuccore = defaultdict(list)
    for protein in genomic_data:
        grouped_by_nuccore[protein['nuccore_id']].append(protein)

    for nuccore_id, proteins in grouped_by_nuccore.items():
        # Sort proteins by their start position
        sorted_proteins = sorted(proteins, key=lambda x: x['start'])

        # Iterate through sorted proteins to find operons
        current_operon = []
        for i in range(len(sorted_proteins)):
            if not current_operon:
                # Start a new operon
                current_operon.append(sorted_proteins[i])
            else:
                last_protein = current_operon[-1]
                # Check if on the same strand and within the threshold
                if (sorted_proteins[i]['strand'] == last_protein['strand'] and
                        (sorted_proteins[i]['start'] - last_protein['stop']) <= threshold):
                    # Add to current operon
                    current_operon.append(sorted_proteins[i])
                else:
                    # Finalize the current operon and assign an operon number
                    operon_number += 1
                    for protein in current_operon:
                        operon_assignments.append({
                            'Protein_ID': protein['Protein_ID'],
                            'nuccore_id': nuccore_id,
                            'start': protein['start'],
                            'stop': protein['stop'],
                            'strand': protein['strand'],  # Include strand information
                            'operon_number': operon_number
                        })
                    # Start a new operon with the current protein
                    current_operon = [sorted_proteins[i]]

        # Finalize the last operon if it exists
        if current_operon:
            operon_number += 1
            for protein in current_operon:
                operon_assignments.append({
                    'Protein_ID': protein['Protein_ID'],
                    'nuccore_id': nuccore_id,
                    'start': protein['start'],
                    'stop': protein['stop'],
                    'strand': protein['strand'],  # Include strand information
                    'operon_number': operon_number
                })

    return operon_assignments

# Write operon results to a file
def write_results(output_file, operon_assignments):
    with open(output_file, 'w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=['Protein_ID', 'nuccore_id', 'start', 'stop', 'strand', 'operon_number'])
        writer.writeheader()
        writer.writerows(operon_assignments)

# Main execution flow
if __name__ == "__main__":
    input_file = 'filtered_file_no_boundaries.csv'  # Replace with your input file
    output_file = 'operon_results.csv'  # Desired output file name
    genomic_data = load_genomic_coordinates(input_file)
    operon_assignments = identify_operons(genomic_data)
    write_results(output_file, operon_assignments)
    print(f"Operon results written to {output_file}")

##Step 3 : combine operon predictions with MMseqs clustering to define operons and representative operons.    
# Load MMseqs clusters from a CSV file and ensure correct headers
def load_clusters(filename):
    clusters = {}
    expected_headers = ['Cluster_ID', 'Protein_ID']  # Define the expected headers

    try:
        # Open the file to read its content
        with open(filename, 'r') as file:
            # Read the first line to check for headers
            first_line = file.readline().strip()
            headers = first_line.split('\t')

            # If the expected headers are not present, prepend them
            if set(expected_headers).issubset(headers):
                # Rewind to read the file as usual
                file.seek(0)
                reader = csv.DictReader(file, delimiter='\t')  # Adjust the delimiter if necessary

                for row in reader:
                    cleaned_row = {k.strip(): v.strip() for k, v in row.items()}  # Clean keys and values
                    if 'Protein_ID' in cleaned_row and 'Cluster_ID' in cleaned_row:
                        clusters[cleaned_row['Protein_ID']] = cleaned_row['Cluster_ID']
            else:
                # Rewind and create a new list of rows with added headers
                file.seek(0)
                new_rows = []
                new_rows.append('\t'.join(expected_headers))  # Add the expected headers

                # Read the rest of the lines and add to new_rows
                for line in file:
                    values = line.strip().split('\t')  # Split by the delimiter
                    if len(values) < len(expected_headers):
                        values.extend([''] * (len(expected_headers) - len(values)))  # Add empty strings for missing fields
                    new_rows.append('\t'.join(values))

                # Write the new rows to a temporary file
                with open('temp_clusters.tsv', 'w') as temp_file:
                    temp_file.write('\n'.join(new_rows))
                
                # Now load the clusters from the temporary file
                with open('temp_clusters.tsv', 'r') as temp_file:
                    reader = csv.DictReader(temp_file, delimiter='\t')
                    for row in reader:
                        cleaned_row = {k.strip(): v.strip() for k, v in row.items()}  # Clean keys and values
                        if 'Protein_ID' in cleaned_row and 'Cluster_ID' in cleaned_row:
                            clusters[cleaned_row['Protein_ID']] = cleaned_row['Cluster_ID']

    except Exception as e:
        print("Error reading the file:", e)

    
    return clusters

# Combine MMseqs data with operon assignments
def combine_operons_with_clusters(operon_assignments, clusters):
    for operon in operon_assignments:
        protein_id = operon['Protein_ID'].strip()  # Clean up whitespace
        cluster_id = clusters.get(protein_id, 'Not Found')  # Get the cluster ID, or 'Not Found' if it doesn't exist
        operon['Cluster_ID'] = cluster_id  # Add the cluster ID to the operon data
  

    return operon_assignments

# Write the combined results to a file, removing rows with 'Not Found' Cluster_ID
def write_combined_results(output_file, operon_assignments):
    valid_assignments = [operon for operon in operon_assignments if operon['Cluster_ID'] != 'Not Found']  # Eliminate rows with 'Not Found'

    # Sort the valid assignments by nuccore_id and start
    valid_assignments.sort(key=lambda x: (x['nuccore_id'], int(x['start'])))

    with open(output_file, 'w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=['Protein_ID', 'nuccore_id', 'start', 'stop', 'strand', 'operon_number', 'Cluster_ID'])
        writer.writeheader()
        writer.writerows(valid_assignments)

# Write representative operons to a separate file
def write_representative_operons(valid_assignments, output_file):
    representative_operons = []
    seen_keys = set()

    for operon in valid_assignments:
        key = (operon['operon_number'], operon['Cluster_ID'])
        if key not in seen_keys:
            representative_operons.append(operon)
            seen_keys.add(key)

    with open(output_file, 'w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=['Protein_ID', 'nuccore_id', 'start', 'stop', 'strand', 'operon_number', 'Cluster_ID'])
        writer.writeheader()
        writer.writerows(representative_operons)
# Write representative operons to a separate file
def write_representative_operons(valid_assignments, output_file):
    representative_operons = []
    seen_keys = set()

    for operon in valid_assignments:
        # Skip operons with "not found" in Cluster_ID
        if operon['Cluster_ID'].lower() == 'not found':  # Check for "not found"
            continue
        
        key = (operon['operon_number'], operon['Cluster_ID'])
        if key not in seen_keys:
            representative_operons.append(operon)
            seen_keys.add(key)

    # Write to output file
    with open(output_file, 'w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=['Protein_ID', 'nuccore_id', 'start', 'stop', 'strand', 'operon_number', 'Cluster_ID'])
        writer.writeheader()
        writer.writerows(representative_operons)

# Main execution flow
if __name__ == "__main__":
    operon_file = 'operon_results.csv'  # Replace with your operon results file
    mmseqs_file = 'MMseqs-clusternumb.tsv'  # Replace with your MMseqs clusters file
    combined_output_file = 'combined_results.csv'  # Desired output file name
    representative_output_file = 'representative_operons.csv'  # File for representative operons

    # Load operon assignments
    operon_assignments = []
    with open(operon_file, 'r') as file:
        reader = csv.DictReader(file)
        operon_assignments = list(reader)

    # Load clusters
    clusters = load_clusters(mmseqs_file)

    # Combine operons with clusters
    combined_results = combine_operons_with_clusters(operon_assignments, clusters)

    # Write combined results to output file
    write_combined_results(combined_output_file, combined_results)

    # Write representative operons to a separate file, excluding 'not found' Cluster_ID
    write_representative_operons(combined_results, representative_output_file)    


