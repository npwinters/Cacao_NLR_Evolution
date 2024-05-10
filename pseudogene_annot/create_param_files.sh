#!/usr/bin/env bash
# Noah Winters
# JAN 2022
# Creates input parameter file for CombinedPseudoWrapper.py
for file in $(ls *_nlr.faa); do
	cat parameter_base |\
	( echo -e "b_out=${file%.faa}.blast\np_seq=${file}\ng_seq=${file%_nlr.faa}.fa"; cat -- ) |\
	( cat --; echo -e "gff=${file%_nlr.faa}_cleaned.gff3\nrepCut=300\nrepDiv=30" ) > ${file%_nlr.faa}_parameter_file.txt

	mkdir ${file%_nlr.faa}
	mv ${file%_nlr.faa}* ${file%_nlr.faa}/
done
