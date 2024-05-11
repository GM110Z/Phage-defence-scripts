# set project directory (update to where your files are)
PROJDIR=$(pwd)

# list the files prefixes to run
FILELIST=$(cat input.txt)

# running PADLOC (assumes padloc is in your PATH)
parallel --j 6 --bar "padloc --cpu 1 --faa $PROJDIR/{}.faa --gff $PROJDIR/{}.gff" ::: $FILELIST
