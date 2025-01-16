**PhiSpy-loop.sh**: runs PhiSpy prophage prediction in loop on multiple files

**OCHRE.py**:OCHRE.py: prOphage defenCe Hotspot boundaRies Extraction: From a precombined database of protein sequences from prophages files, script gets all proteins specified by user in sys.argv[1] and retrieves all info from NCBI. If user want, they can to specify as sys.argv[2] a file containing a list of Assembly IDs for the script to keep when dropping duplicates. Script produces output to run in https://github.com/GCA-VH-lab/FlaGs2

**PARSEC**:Finds prophages in genomes, groups them with fastANI It needs PhySpy, and faidx, fastANI and ANIClustermap(https://github.com/moshi4/ANIclustermap/tree/main)

**fastani-to-clusters.py**: runs as fastani-to-clusters.oy <threshold for clustering (float)> Script from https://github.com/moshi4/ANIclustermap/tree/main), needs pandas
Run it like: python fastani-to-clusters.py <anicluster-matrix> 95 (threshold for similarity) <output-file-name>


**ViralTreeAnno** : R code that uses ggtree to annotate a VIPTree or any tree with user input data in table format

**find-defence-prophages.py** compares output of fastANI-clusters from all prophages vs prophages that had defences, to select clusters from the table of All prophages that have members with a defence system

**probeV2.1.py**: PRotein cOlocalisation By hmmEr: Uses HMM models of an upstream and downstream protein to extract genomic regions. It allows an intergenic distance of 10kb.
sys.argv[1]hmm model of boundary 1
sys.argv[2]:hmm model of boundary 2
sys.argv[3]:directory with faa and gb files 


**PECAN.py**: ProtEin Colocalization by ANnotations: Uses the annotations in genbank 'product' to find desired proteins. It finds only those that colocalise <10000kb distance:
sys.argv[1]:Annotation search term1
sys.argv[2]:Annotation search term2
sys.argv[3]:Allowed intergenic distance between search terms
sys.argv[4]:path/to/genbank/files
sys.argv[5]:Name of outputfile

**PHORIFIC.py** It needs local genbank files (*.gb) to build a table of coordinates for each proteins. Then based on nuccore ID and coordinates, predicts which proteins are in operons. COmbines these data with MMSeqs clustering (run separately). This allows to find operons that occurr more than once and pick one representative each for further analysis 
NOTE:FOR NOW EDITING OF INPUT FILES  NAMES AND PATH AND OPERON DISTANCE NEEDS TO BE MANUALLY CHANGED --FIX WILL FOLLOW SOON WHEN I HAVE DECIDED IF THE SCRIPT IS COMPLETE AND GOOD LIKE THIS

**edison.py**: Edison (hiddEn moDel proteIn familieS identificatiON)Runs hmmscan vs a PFAM database and processes data to only report hits with e-value <0.01 

**SantasHelper** Uses pandas to combine outputs of PHORIFIC, PADLOC and EDISON. Remember to run padloc separately.


**ICARUS.py** Uses python to download representative sequences for each predicted operons efetch and then runs cblaster to determine if this is effectively an operon (occurs in multiple genomes). Cblaster run is times so it does not overwhelm blast. Originally written for the project https://github.com/GM110Z/Shield-subtypes, this version of ICARUS leaves out the DUF2130 requirement by commenting it out. 

**get_binary_files.sh** extracts only the binary files produced by cblaster if they have 100 or more lines (thus 100 or more operons) and moves them in a different folder
