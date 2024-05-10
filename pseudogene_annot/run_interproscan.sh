#!/usr/bin/env bash
# Noah Winters
# JAN 2022
# Runs interproscan 
for file in $(ls *.faa); do
        bash interproscan.sh -i $file -f tsv -dp -appl Pfam -o ${file%.faa}_pfam.tsv -T /scratch/users/noah/software/PseudogenePipeline/genotype_files/temp
done
