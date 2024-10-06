import os
from Bio import SeqIO
import subprocess

# Run PFAM annotation using HMMER
def annotate_with_pfam(fasta_file, pfam_db, output_file):
    """Annotate protein sequences with PFAM domains using HMMER."""
    if not os.path.exists(fasta_file):
        raise FileNotFoundError(f"FASTA file '{fasta_file}' not found.")
    if not os.path.exists(pfam_db):
        raise FileNotFoundError(f"PFAM database file '{pfam_db}' not found.")
    
    command = f"hmmscan --domtblout {output_file} {pfam_db} {fasta_file}"
    try:
        subprocess.run(command, shell=True, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running hmmscan: {e.stderr}")
        raise

# Parse the HMMER output file
def parse_hmmer_output(output_file):
    """Parse HMMER domain table output to extract Pfam domains, keeping only the best E-value per target name."""
    pfam_domains = {}
    with open(output_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split()
            if len(fields) < 23:
                print(f"Skipping malformed line: {line}")
                continue
            
            query_name = fields[0]  # Protein ID (query name)
            target_name = fields[3]  # Pfam domain (target name)
            accession = fields[1]     # Pfam accession
            e_value = float(fields[12])  # E-value
            description = ' '.join(fields[22:]) if len(fields) > 22 else 'Unknown'

            # Only consider entries with E-value <= 0.01
            if e_value <= 0.01:
                # Keep only the best (lowest) E-value for each unique target name
                if target_name not in pfam_domains or e_value < pfam_domains[target_name][2]:
                    pfam_domains[target_name] = (query_name, accession, e_value, description)

    return pfam_domains

# Save the parsed Pfam domains to a TSV file
def save_to_tsv(pfam_domains, output_file):
    """Save the parsed Pfam domains and descriptions to a TSV file, ensuring no duplicates for the target name."""
    with open(output_file, 'w') as f:
        f.write("target name\taccession\tquery name\tE-value\tdescription of target\n")
        for target, values in pfam_domains.items():
            query_name, accession, e_value, description = values
            f.write(f"{target}\t{accession}\t{query_name}\t{e_value:.2E}\t{description}\n")

# Example usage
fasta_file = "database.fa"
pfam_db = "Pfam-A.hmm"
output_file = 'pfam_annotations.txt'

# Step 1: Run hmmscan to create the output file
annotate_with_pfam(fasta_file, pfam_db, output_file)

# Step 2: Parse the HMMER output file (only if the previous step succeeds)
parsed_domains = parse_hmmer_output(output_file)

# Step 3: Save the parsed domains to a new TSV file
save_to_tsv(parsed_domains, 'filtered_pfam_domains_with_desc.tsv')

