import sys
from Bio import SeqIO, SeqFeature
import csv
import glob
import os
import pandas as pd
import numpy as np


print('Working on the genomes...')
#combine gbk file in working dir
for f in glob.glob("*.gb"):
         os.system("cat "+f+" >> OutFile.gb")

# Parse combined records
combined_files = SeqIO.parse("OutFile.gb", format="genbank")

with open('prefiltered.txt', mode='w') as parsed_output:   #write table for Suppl.data with hitsID, scores and lenght
                parsed_output = csv.writer(parsed_output, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)
                parsed_output.writerow(['Protein_ID','product', 'molecule', 'start', 'end', 'strand'])
# Print the info of each record
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
df5.to_csv("gggenes-input.txt", "\t", index=False, header=True)

table1 = pd.read_csv('gggenes-input.txt', sep="\t")
table2 = pd.read_csv('Knonw-and-New-candidates-in-hotspot.csv', sep=",")

# Merge the two tables on 'Protein_ID', keeping all rows from table1 (left join)
merged_df = pd.merge(
    table1,
    table2[['Protein_ID', 'Description_of_Target', 'system']],
    on='Protein_ID',
    how='left'
)

# Save the merged DataFrame to a new CSV file
merged_df.to_csv('merged_table.csv', index=False)
                
#removing combined file                
os.system("rm OutFile.gb") 
os.system("rm prefiltered.txt") 
os.system("rm gggenes-input.txt") 
