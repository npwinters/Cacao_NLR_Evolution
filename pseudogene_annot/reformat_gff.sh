#!/usr/bin/env bash
# Noah Winters
# JAN 2022
# Use to add "longest" and "Parent" to GFF3s so that they will run with PseudogenePipeline
# This assumes your GFF3 contains only the longest isoform for each gene (usually .1)

awk -F "\t" -v OFS="\t" '{
if($3=="mRNA"){gsub("Parent","longest=1;Parent",$9); split($9,ID,";"); split(ID[1],pac,"="); $9=$9";pacid="pac[2]};
if($3=="CDS"||$3=="exon"||$3~/UTR$/){split($9,ID,";"); split(ID[2],pac,"="); $9=$9";pacid="pac[2]}
}1' $1
