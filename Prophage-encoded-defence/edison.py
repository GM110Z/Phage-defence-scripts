#run PFAM
import os
from Bio import SeqIO
import subprocess

def annotate_with_pfam(fasta_file, pfam_db, output_file):
    """Annotate protein sequences with PFAM domains using HMMER."""
    command = f"hmmscan --domtblout {output_file} {pfam_db} {fasta_file}"
    subprocess.run(command, shell=True, check=True)

def parse_pfam_output(pfam_output_file):
    """Parse HMMER output file to extract PFAM domains."""
    annotations = {}
    with open(pfam_output_file, 'r') as f:
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue
            parts = line.split()
            seq_id = parts[0]
            pfam_id = parts[3]
            if seq_id not in annotations:
                annotations[seq_id] = []
            annotations[seq_id].append(pfam_id)
    return annotations

# Example usage
fasta_file = "database.faa"  # Your input FASTA file
pfam_db = "Pfam-A.hmm"  # Your PFAM HMM database
pfam_output_file = "pfam_annotations.txt"  # Output file for PFAM annotations

annotate_with_pfam(fasta_file, pfam_db, pfam_output_file)
annotations = parse_pfam_output(pfam_output_file)

##process the output file
def parse_hmmer_output(output_file):
    """Parse HMMER domain table output to extract Pfam domains and sequence descriptions, keeping only the best E-value per sequence."""
    pfam_domains = {}
    sequence_descriptions = {}
    with open(output_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split()
            if len(fields) < 22 or fields[0].startswith('#'):
                continue
            sequence_id = fields[0]
            domain_id = fields[3]
            domain_name = fields[4]
            e_value = float(fields[12])
            description = ' '.join(fields[22:]) if len(fields) > 22 else 'Unknown'
            
            if e_value <= 0.01:  # Exclude domains with E-value > 0.01
                # If this sequence_id is new, or this e_value is better than the existing one, keep it
                if (sequence_id not in pfam_domains) or (e_value < pfam_domains[sequence_id][2]):
                    pfam_domains[sequence_id] = (domain_id, domain_name, e_value)
                    sequence_descriptions[sequence_id] = description
    return pfam_domains, sequence_descriptions

def save_to_tsv(pfam_domains, sequence_descriptions, output_file):
    """Save parsed Pfam domains and sequence descriptions to a TSV file."""
    with open(output_file, 'w') as f:
        f.write("Sequence_ID\tDescription_of_Target\tDomain_ID\tDomain_Name\tE_Value\n")
        for seq_id, domain in pfam_domains.items():
            description = sequence_descriptions.get(seq_id, 'Unknown')
            f.write(f"{seq_id}\t{description}\t{domain[0]}\t{domain[1]}\t{domain[2]}\n")

# Example usage:
output_file = 'filtered_pfam_domains_with_desc.tsv'
parsed_domains, sequence_desc = parse_hmmer_output('pfam_annotations.txt')
save_to_tsv(parsed_domains, sequence_desc, output_file)
print(f"Filtered Pfam domains with descriptions saved to {output_file}")
