#!/usr/bin/env bash
# Noah Winters
# JAN 2022
# Extract chromosome ID, TE ID, and TE type from RepeatMasker .out file
for file in $(ls *Theobroma_cacao*.out); do 
	cat $file | tail +4 | tr -s " " | sed 's/^ //g' | cut -d ' ' -f 10,11 | sed 's/ /	/g' | sort -k1 | uniq > ${file}.type.txt 
done

# Make tmp gff that has position of each feature 
for file in $(ls *gff); do 
	cat $file | cut -f 1,4,5,7,9 | tail +4 | sed 's/Target "Motif\://g' | sed 's/" .*//g' | awk '{print $5"\t"$1"\t"$2"\t"$3"\t"$4}' > ${file%.gff}.tmp.gff
done


# Merge TE description with gff position information
for file in $(ls *.type.txt); do 
	bash innerJoin_stdOut.sh $file ${file%.type.txt}.tmp.gff | tail +2 | cut -f 4,5,6,7,1,2 | awk '{print $3"\t"$4"\t"$5"\t"$1"\t"$2"\t"$6}' > ${file}.TE.pos.txt
done

