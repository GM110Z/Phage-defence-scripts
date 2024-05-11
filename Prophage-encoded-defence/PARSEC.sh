
#!/bin/bash


VAR=$(echo *.gbff | xargs ls)

for file in ${VAR}

        do
                var2=$(echo $file)
                PhiSpy.py ${file} -p ${var2}  -o ${var2}_output_folder --output_choice 31 --threads 8;


        done
#create fasta database
cwd=$(pwd)
find . -name '*.gbff_phage.fasta' -exec cp {} ${cwd}/ \

#Split multifasta of multiple prophages

VAR2=$(ls *.gbff_phage.fasta)

for i in ${VAR2}
   do
      faidx --split-files ${i}
   done

rm *.fai
rm *.gbff_phage.fasta


#cluster prophages by similarity with fastANI and ANIClustermap
ANIclustermap -i $1 -o $2 --cmap_colors white,orange,red --fig_height 200 --fig_width 50


#if you want to run normal fastANI uncomment below
#ls *.fasta > genome-paths.txt
#fastANI --ql genome-paths.txt --rl genome-paths.txt -t 12 -o fastANI-out-raw.tsv
#awk '$1 != $2 ' fastANI-out-raw.tsv | awk ' { palign = $4 / $5 * 100 } { print $0"\t"palign } ' | sort -nr -k 6,6 -k 4,4 > fastANI-out.tsv
#cat <( printf "query\treference\tANI\tnum_frags_mapped\ttotal_query_frags\tpercent_aligned\n" ) fastANI-out.tsv > fastANI-out-with-header.tsv


