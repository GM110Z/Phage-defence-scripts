import os
import sys
import glob
import csv
from Bio import SeqIO

# Function to search for specific product descriptions and write to CSV
def search_genbank_products(genbank_directory, output_csv):
    # Search terms
    search_terms = [sys.argv[1], sys.argv[2]]

    # Open CSV for writing
    with open(output_csv, mode='w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(['ID', 'Start', 'Stop'])  # CSV header

        # Process each GenBank file
        for genbank_file in glob.glob(f"{genbank_directory}/*.gbff"):
            print(f"Processing GenBank file: {genbank_file}")
            
            # Use SeqIO.parse to handle multiple records in a file
            for genbank_record in SeqIO.parse(genbank_file, "genbank"):

                # Collect protein locations
                protein_locations = []
                for feature in genbank_record.features:
                    if feature.type == "CDS" and "product" in feature.qualifiers:
                        product_desc = feature.qualifiers["product"][0].lower()
                        if any(term in product_desc for term in search_terms):
                            location = feature.location
                            protein_locations.append((product_desc, location.start, location.end, genbank_record.id))

                # Write to CSV only if proteins are close to each other
                for i, (product1, start1, end1, chrom_id1) in enumerate(protein_locations):
                    for j, (product2, start2, end2, chrom_id2) in enumerate(protein_locations):
                        if i != j and chrom_id1 == chrom_id2:
                            if abs(start1 - end2) <= sys.argv[3] or abs(start2 - end1) <= sys.argv[3]:
                                combined_start = min(start1, start2)
                                combined_end = max(end1, end2)
                                writer.writerow([chrom_id1, combined_start, combined_end])
                                print(f"Written coordinates: {chrom_id1}, {combined_start}-{combined_end}")

    print(f"Results written to {output_csv}")

# Example call to search the GenBank files
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python search_genbank_products.py <genbank_directory> <output_csv>")
        sys.exit(1)

    genbank_directory = sys.argv[4]
    output_csv = sys.argv[5]
    search_genbank_products(genbank_directory, output_csv)

