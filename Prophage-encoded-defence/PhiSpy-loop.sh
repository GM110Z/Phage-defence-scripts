#!/bin/bash


VAR=$(echo *.gbff | xargs ls)

for file in ${VAR}

        do
                var2=$(echo $file)
                PhiSpy.py ${file} -p ${var2}  -o ${var2}_output_folder --threads 8 --output_choice 31;


        done
