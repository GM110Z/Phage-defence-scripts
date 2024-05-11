#!/usr/bin/env python
import pandas as pd
import numpy
import Bio
import sys
from Bio import SeqIO

for record in SeqIO.parse('database.fa', 'fasta'):
    if sys.argv[1] in record.description:
        print (record.id,file=open('efetch-input.txt','a'))


list_of_accession = []
with open ('efetch-input.txt', 'r', encoding='utf-8-sig') as csvfile:
    efetchin=csv.reader(csvfile, delimiter = ',')
    for row in efetchin:
        list_of_accession.append(str(row[0]))
        
outfile=open('efetch_output.txt', mode = 'wb')
input_handle = Entrez.efetch(db="protein", id= list_of_accession, rettype="ipg", retmode="tsv")
for line in input_handle:
    outfile.write(line)
input_handle.close()
outfile.close()


df = pd.read_csv("efetch_output.txt",sep="\t", low_memory=False)
df.columns=['ID','Source', 'Nucleotide Accession', 'Start', 'Stop', 'Strand', 'Protein', 'Protein Name', 'Organism', ' Strain', 'Assembly']
df2 = pd.read_csv(sys.argv[2],sep="\t", low_memory=False)
lists=df2.to_numpy().flatten().tolist()
df3=df.loc[df['Assembly'].isin(lists)]
df4 = df3.drop_duplicates(subset=['Protein','Assembly'],keep='first')
df4.to_csv('selected-proteins-coord.tsv', '\t', header=True)
df4.sort_values(by = "Organism", axis=0, ascending=True, inplace=False).to_csv('flags-input.tsv', '\t', header=False, index= False, columns = ['Assembly', 'Protein'])
