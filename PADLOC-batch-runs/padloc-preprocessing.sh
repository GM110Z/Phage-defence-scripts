#!/bin/bash
#convert gb to gff3
bp_genbank2gff3 *.gb

echo 'cleaning gff files'
#clean gff for padloc
VAR4=$(echo *.gff | xargs ls)

for file in ${VAR4}

        do
           sed '/^##FASTA/Q' ${file} >${file}_noseq.gff
        done

#conver gb to protein fasta
VAR4=$(echo *.gb | xargs ls)

for f in ${VAR4}
   do
      python gbktofaa.py  ${f}
   done
