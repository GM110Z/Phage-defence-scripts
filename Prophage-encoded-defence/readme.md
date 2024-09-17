**PhiSpy-loop.sh**: runs PhiSpy prophage prediction in loop on multiple files

**OCHRE.py**:OCHRE.py: prOphage defenCe Hotspot boundaRies Extraction: From a precombined database of protein sequences from prophages files, script gets all proteins specified by user in sys.argv[1] and retrieves all info from NCBI. If user want, they can to specify as sys.argv[2] a file containing a list of Assembly IDs for the script to keep when dropping duplicates. Script produces output to run in https://github.com/GCA-VH-lab/FlaGs2

**PARSEC**:Finds prophages in genomes, groups them with fastANI It needs PhySpy, and faidx, fastANI and ANIClustermap(https://github.com/moshi4/ANIclustermap/tree/main)

**fastani-to-clusters.py**: runs as fastani-to-clusters.oy <threshold for clustering (float)> Script from https://github.com/moshi4/ANIclustermap/tree/main), needs pandas


**ViralTreeAnno** : R code that uses ggtree to annotate a VIPTree or any tree with user input data in table format

**find-defence-prophages.py** compares output of fastANI-clusters from all prophages vs prophages that had defences, to select clusters from the table of All prophages that have members with a defence system

**probeV2.1.py**: PRotein cOlocalisation By hmmEr: Uses HMM models of an upstream and downstream protein to extract genomic regions. 1st arg: directory with .faa files 2nd arg : directory to genbank files. Runs without wrapper. 

**pecan.py**"ProtEin Colocalization by ANnotations: Uses the annotations in genbank 'product' to find desired proteins. It finds only those that colocalise <10000kb distance:
sys.argv[1]:Annotation search term1
sys.argv[2]:Annotation search term2
sys.argv[3]:Allowed intergenic distance between search terms
sys.argv[4]:path/to/genbank/files
sys.argv[5]:Name of outputfile
