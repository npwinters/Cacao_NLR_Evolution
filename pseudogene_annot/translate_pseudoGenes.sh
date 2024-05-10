#!/usr/bin/env bash
# Noah Winters
# JAN 2022
# Translates annotated pseudogenes in all 6 frames
FILE=$1
cat $FILE | sed '/^#/d' | translate6frames.sh -in=stdin.fa  out=stdout.fa | sed 's/\*/X/g' | sed 's/\./X/g' | sed 's/ fr/_fr/g' > ${FILE}.faa


