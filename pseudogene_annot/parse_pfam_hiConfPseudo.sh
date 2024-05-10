#!/usr/bin/env bash
# Noah Winters
# JAN 2022
# Parses high confidence pseudogene identifications from Pfam
FILE=$1
cat $FILE | grep -f nlr_pfam_id.txt | cut -f 1 | sed 's/_fr.*//g' | sort | uniq > ${FILE%_pfam.tsv}_cleaned_pfam.txt
