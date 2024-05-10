#!/usr/bin/env bash
# Noah Winters
# JAN 2022
# Run tblastn
for file in $(ls *_nlr.faa); do
	tblastn -query $file -db blastdb/${file%_nlr.faa}.fa -evalue 1e-20 -outfmt 6 -out ${file%.faa}.blast
done
