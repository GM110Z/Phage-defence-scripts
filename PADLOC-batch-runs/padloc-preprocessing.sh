#!/bin/bash
#convert gb to gff3
VAR4=$(echo *.gb | xargs ls)

for f in ${VAR4}
   do
      python gbktogff.py  ${f}
   done


#conver gb to protein fasta
VAR4=$(echo *.gb | xargs ls)

for f in ${VAR4}
   do
      python gbktofaa.py  ${f}
   done
