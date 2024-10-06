import os
import sys
import glob
import subprocess as sp
import csv
import time
from Bio import SearchIO
from Bio import SeqIO

# HMMER search function with adjustable thresholds
def run_hmmer(model, fasta_file, output_file, threshold=10):
    print(f"Running HMMER on {fasta_file} with model {model}...")
    hmmsearch_cmd = f"hmmsearch -T {threshold} --incT {threshold} -o log --domtblout {output_file} {model} {fasta_file}"
    sp.run(hmmsearch_cmd, shell=True)

# Function to parse HMMER output and write to CSV
def parse_hmmer_output(hmm_output_file, csv_output_file):
    if not os.path.exists(hmm_output_file):
        print(f"Warning: {hmm_output_file} does not exist.")
        return
    with open(hmm_output_file, newline='') as input_file:
        for qresult in SearchIO.parse(input_file, 'hmmscan3-domtab'):
            hits = qresult.hits
            if hits:
                with open(csv_output_file, mode='w', newline='') as parsed_output:
                    writer = csv.writer(parsed_output, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
                    for hit in hits:
                        writer.writerow([hit.id, hit.description, hit.evalue, hit.bitscore, hit.seq_len])
                print(f"Hits found and written to {csv_output_file}")
            else:
                print(f"No hits found in {hmm_output_file}")

# Function to get protein location from GenBank
def get_protein_location(genbank_file, protein_id):
    protein_location = None
    # Parse the GenBank file to find the record
    for genbank_record in SeqIO.parse(genbank_file, "genbank"):
        for feature in genbank_record.features:
            if feature.type == "CDS" and protein_id in feature.qualifiers.get("protein_id", []):
                protein_location = (feature.location.start, feature.location.end, genbank_record.id)
                return protein_location  # Return as soon as we find a match
    return None, None, None

# Function to extract protein IDs from CSV files
def extract_protein_id(csv_file):
    if not os.path.exists(csv_file):
        print(f"Warning: {csv_file} does not exist.")
        return None
    with open(csv_file, 'r', encoding='utf-8-sig') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            return row[0]  # Assuming the first column has the protein ID

# Main processing function
def process_files(model1, model2, protein_directory, genbank_directory):
    fasta_files = glob.glob(f"{protein_directory}/*.faa")
    genbank_files = glob.glob(f"{genbank_directory}/*.gbff")

    # Dictionary for GenBank file lookup
    genbank_dict = {os.path.splitext(os.path.basename(f))[0]: f for f in genbank_files}

    print(f"Found {len(fasta_files)} FASTA files and {len(genbank_files)} GenBank files.")

    for fasta_file in fasta_files:
        base_name = os.path.splitext(os.path.basename(fasta_file))[0]

        # Check if matching GenBank file exists
        if base_name in genbank_dict:
            genbank_file = genbank_dict[base_name]
            print(f"Processing: {fasta_file} and {genbank_file}")

            # Run HMMER for both models
            run_hmmer(model1, fasta_file, "HMM_output-search1.txt", threshold=10)
            run_hmmer(model2, fasta_file, "HMM_output-search2.txt", threshold=10)

            # Parse the HMMER output files
            parse_hmmer_output("HMM_output-search1.txt", "protein1.csv")
            parse_hmmer_output("HMM_output-search2.txt", "protein2.csv")

            # Extract protein IDs from CSVs
            protein_id1 = extract_protein_id("protein1.csv")
            protein_id2 = extract_protein_id("protein2.csv")

            if protein_id1 and protein_id2:
                # Get protein locations from GenBank
                start1, end1, chrom_id1 = get_protein_location(genbank_file, protein_id1)
                start2, end2, chrom_id2 = get_protein_location(genbank_file, protein_id2)

                # Combine coordinates if proteins are on the same chromosome and within 10,000 base pairs
                if all([start1, end1, chrom_id1, start2, end2, chrom_id2]) and chrom_id1 == chrom_id2:
                    if abs(start1 - end2) <= 10000 or abs(start2 - end1) <= 10000:
                        combined_start = min(start1, start2)
                        combined_end = max(end1, end2)
                        with open("Brex-coordinates.txt", "a") as f:
                            f.write(f"{chrom_id1}\t{combined_start}\t{combined_end}\n")
                        print(f"Written coordinates: {chrom_id1}, {combined_start}-{combined_end}")
                    else:
                        print(f"Proteins are on the same chromosome but too far apart: {protein_id1}, {protein_id2}")
                else:
                    print(f"Proteins are on different chromosomes or not found: {protein_id1}, {protein_id2}")
            else:
                print("One or both protein IDs could not be extracted.")

            # Clean up temporary files
            if os.path.exists("HMM_output-search1.txt"):
                os.remove("HMM_output-search1.txt")
            if os.path.exists("HMM_output-search2.txt"):
                os.remove("HMM_output-search2.txt")
            if os.path.exists("protein1.csv"):
                os.remove("protein1.csv")
            if os.path.exists("protein2.csv"):
                os.remove("protein2.csv")
            time.sleep(2.5)
        else:
            print(f"No matching GenBank file found for {fasta_file}.")


# Example call to process the files
if __name__ == "__main__":
    model1 = sys.argv[1]
    model2 = sys.argv[2]
    protein_directory = sys.argv[3]
    genbank_directory = sys.argv[3]
    process_files(model1, model2, protein_directory, genbank_directory)
