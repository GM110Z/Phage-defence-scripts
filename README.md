**probe.py**: PRotein cOlocalisation By hmmEr: Uses HMM models of an upstream and downstream protein to extract genomic regions.
1st arg: directory with .faa files
2nd arg : genbank file

Run the script within a bash wrapper:
**wrapper.sh** $1
Where $1: a text file listing all genomes IDs (without extension)
