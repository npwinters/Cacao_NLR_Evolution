#!/usr/bin/env bash
# Noah Winters
# JAN 2022
# Runs PseudogenePipeline from the Shiu lab
while read genotype; do
	echo $genotype
	python /scratch/users/noah/software/PseudogenePipeline/_wrapper_scripts/CombinedPseudoWrapper.py ${genotype}_parameter_file.txt

	mv _results ${genotype}_results
	mv _intermediate ${genotype}_intermediate
	mv _logs ${genotype}_logs

done < $1

rm *gff3.gene4col *.noAlt GFFFeatures.out
