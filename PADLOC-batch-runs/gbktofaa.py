#!/usr/bin/env python
# coding: utf-8
from Bio import SeqIO
import sys
file_name = sys.argv[1]

# stores all the CDS entries
all_entries = []

with open(file_name, 'r') as GBFile:

    GBcds = SeqIO.InsdcIO.GenBankCdsFeatureIterator(GBFile)

    for cds in GBcds:
        if cds.seq is not None:
            cds.id = cds.id     #can also just use cds.id
            cds.description = ''
            all_entries.append(cds)


# write file
SeqIO.write(all_entries, '{}.fasta'.format(file_name[:-3]), 'fasta')
