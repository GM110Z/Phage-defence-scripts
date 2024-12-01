#!/usr/bin/env python
from Bio import SeqIO
import sys

def gbk_to_gff3(gbk_file, output_file):
    # Open the GenBank file
    with open(gbk_file, "r") as handle, open(output_file, "w") as out_handle:
        # Write GFF header (optional)
        out_handle.write("##gff-version 3\n")
        for record in SeqIO.parse(handle, "genbank"):
            out_handle.write(f"##sequence-region {record.id} 1 {len(record)}\n")

            # Iterate over features and convert them to GFF3 format
            for feature in record.features:
                if feature.type not in ['source', 'repeat_region', 'low_complexity', 'rRNA', 'tRNA', 'misc_RNA']:
                    # Default GFF3 format structure
                    seq_id = record.id
                    source = "genbank"
                    feature_type = feature.type
                    start = feature.location.start + 1  # GFF is 1-based
                    end = feature.location.end
                    score = "."
                    strand = "+" if feature.strand == 1 else "-" if feature.strand == -1 else "."
                    phase = "."  # Phase is only relevant for CDS features
                    attributes = []

                    # Add additional attributes
                    if "gene" in feature.qualifiers:
                        attributes.append(f"ID={feature.qualifiers['gene'][0]}")
                    if "product" in feature.qualifiers:
                        attributes.append(f"Name={feature.qualifiers['product'][0]}")
                    if "locus_tag" in feature.qualifiers:
                        attributes.append(f"locus_tag={feature.qualifiers['locus_tag'][0]}")

                    # Join attributes with semicolons
                    attributes_str = ";".join(attributes) if attributes else "."

                    # Format the feature line in GFF3 format
                    gff_line = "\t".join([seq_id, source, feature_type, str(start), str(end), score, strand, phase, attributes_str])
                    out_handle.write(gff_line + "\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python gbk_to_gff3.py <genbank_file> <output_gff_file>")
        sys.exit(1)
    
    gbk_file = sys.argv[1]
    output_file = sys.argv[2]
    gbk_to_gff3(gbk_file, output_file)
