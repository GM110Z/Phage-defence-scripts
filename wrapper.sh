FILELIST=$1

for i in ${FILELIST}
   do
      python probe.py ${i}.faa ${i}.gb
   done
