#!/bin/bash 
echo "Downloading and measuring genomes"
#download genomes where hits were found
parallel -a fadix_file.tsv -C '\t' -j8 --delay 1.5 "efetch -db nuccore -id {1} -seq_start {2} -seq_stop {3} -format fasta >{1}.fasta"  


#get sizes of each chromosome with faidx
VAR2=$(echo *.fasta | xargs ls)
for y in ${VAR2}
   do
      faidx ${y} -i chromsizes > ${y}.size
   done
#remove whole chromosomes
rm *.fasta
rm *.fai

cat *.size>>chromosome.sizes

echo "Getting neighbourhood"

#get flanking regions
##split the bed files group in single files
awk '{print $1 "\t" $2 "\t" $3 > $1 ".bed"}' Bed_file.tsv

##run bedslop getting 25 bp  each side
VAR3=$(echo *.bed | xargs ls)
for i in ${VAR3}
   do
      bedtools slop -i ${i} -g chromosome.sizes -b 50000 > ${i}.slop
   done
cat *.slop >> regions.boundaries
rm *.bed
rm *.slop
rm *.size
rm chromosome.sizes

#download genbank files of selected regions
parallel -a regions.boundaries -C '\t' -j8 --delay 0.4 "efetch -db nuccore -id {1} -seq_start {2} -seq_stop {3} -format gbwithparts >{1}.gb"  

echo "Converting gbk to fasta protein"
#convert gbk to fasta protein
VAR4=$(echo *.gb | xargs ls)

for f in ${VAR4}
   do
      python gbktofaa.py  ${f}
   done


#remove empty sequence files 
find . -size 0 -exec rm {} \;
#adds underscore to spaces-only works for mac
sed -i '' 's/ /_/g' *.fasta


cat *.fasta >> database.fa

echo "Clustering protein families"
#cluster proteins by similarity
mmseqs easy-linclust database.fa clusterRes tmp --cov-mode 0 -c 0.5 --min-seq-id 0.6 --threads 8
