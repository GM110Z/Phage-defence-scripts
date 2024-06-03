#HMMER run
import subprocess as sp
import csv
import sys
import os
from Bio import SearchIO
from Bio import SeqIO
import time


#looking for brex proteins in each genomes
print('running hmmer...')

hmmsearch = "hmmsearch" + " " + "-T 30" + " "+ "--incT 30" +" "+ "-o log"+ " --domtblout"

if __name__ == '__main__':
    model1= 'pglx.hmm'
    model2= 'pglz.hmm'
    protein_directory =sys.argv[1]

    # Make search for protein1
    hmmsearch_cmd = f"{hmmsearch} HMM_output-search1.txt {model1} {protein_directory}"
    sp.run(hmmsearch_cmd, shell=True)

    with open('HMM_output-search1.txt', newline='') as input:
        for qresult in SearchIO.parse(input, 'hmmscan3-domtab'):
            query_id = qresult.id
            query_len = qresult.seq_len
            hits = qresult.hits
            num_hits = len(hits)
            hits_len = qresult.seq_len

            if num_hits > 0:
                with open('protein1.csv', mode='w', newline='') as parsed_output:
                    parsed_output = csv.writer(parsed_output, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
                    for i in range(0, num_hits):
                        hit_evalue = hits[i].evalue
                        hit_bit_score = hits[i].bitscore
                        hit_accession = hits[i].id
                        hit_length = hits[i].seq_len
                        hit_description = hits[i].description
                        parsed_output.writerow([hit_accession, hit_description, hit_evalue, hit_bit_score, hit_length])

    # Make search for protein2
    hmmsearch_cmd = f"{hmmsearch} HMM_output-search2.txt {model2} {protein_directory}"
    sp.run(hmmsearch_cmd, shell=True)

    with open('HMM_output-search2.txt', newline='') as input:
        for qresult in SearchIO.parse(input, 'hmmscan3-domtab'):
            query_id = qresult.id
            query_len = qresult.seq_len
            hits = qresult.hits
            num_hits = len(hits)
            hits_len = qresult.seq_len

            if num_hits > 0:
                with open('protein2.csv', mode='a', newline='') as parsed_output:
                    parsed_output = csv.writer(parsed_output, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
                    for i in range(0, num_hits):
                        hit_evalue = hits[i].evalue
                        hit_bit_score = hits[i].bitscore
                        hit_accession = hits[i].id
                        hit_length = hits[i].seq_len
                        hit_description = hits[i].description
                        parsed_output.writerow([hit_accession, hit_description, hit_evalue, hit_bit_score, hit_length])
#process HMMER runs and get the coordinates
print('getting those locations...')

def get_protein_location(genbank_file, protein_id):
    # Parse the GenBank file
    genbank_record = SeqIO.read(genbank_file, "genbank")

    # Iterate through the features in the record
    for feature in genbank_record.features:
        if feature.type == "CDS":  # We're interested in protein-coding features
            # Check if the protein ID matches
            if protein_id in feature.qualifiers.get("protein_id", []):
                location = feature.location
                strand = location.strand
                if strand == 1:
                    start, end = location.start, location.end
                elif strand == -1:
                    start, end = location.end, location.start
                else:
                    raise ValueError(f"Invalid strand value: {strand}")

                chrom_id = genbank_record.id
                return start, end, chrom_id

    # If no matching protein ID is found, return None
    return None, None, None

#and now get the brex only coordinate
genbank_file = sys.argv[2]  #this will need to be a sys arg2
protein_id = ''
with open ('protein1.csv', 'r', encoding='utf-8-sig') as csvfile:
    efetchin=csv.reader(csvfile, delimiter = ',')
    for row in efetchin:
        protein_id += row[0]


protein_id2 =''
with open ('protein2.csv', 'r', encoding='utf-8-sig') as csvfile:
    efetchin=csv.reader(csvfile, delimiter = ',')
    for row in efetchin:
        protein_id2 += row[0] + ''
start1, end1, chrom_id1 = get_protein_location(genbank_file, protein_id)
start2, end2, chrom_id2 = get_protein_location(genbank_file, protein_id2)

if start1 is not None and end1 is not None and chrom_id1 is not None and start2 is not None and end2 is not None and chrom_id2 is not None:
    if chrom_id1 == chrom_id2:
        combined_start = min(start1, start2)
        combined_end = max(end1, end2)
        f = open("Brex-coordinates.txt", "a")
        print(f"{chrom_id1}\t{combined_start}\t{combined_end}", file=f)
    else:
        print(f"Proteins {protein_id} and {protein_id2} are on different chromosomes {chrom_id1} and {chrom_id2}.")
else:
    print("One or both of the proteins were not found in the GenBank file.")

time.sleep(2.5)
#now remove the hmmfiles for a new run 
os.system("rm HMM_output-search1.txt") 
os.system("rm HMM_output-search2.txt")
os.system("rm protein1.csv")
os.system("rm protein2.csv")
