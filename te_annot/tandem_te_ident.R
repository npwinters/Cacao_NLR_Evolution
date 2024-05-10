#!/usr/bin/env bash
# Noah Winters
# JAN 2022
# Identifiess TE arrays across a genome
library(tidyverse)
library(intervals)
library(IRanges)
library(GenomicRanges)
library(plyranges)

setwd("/Users/noah/GradSchool/Research/Projects/NSF_PlantGenome/Comparative_Genomics/ideogram_files/SCA6/")
gff_te <- read.table("Theobroma_cacao_SCA-6_chr.v1.0.fa.out.TE.gff", header = FALSE, sep = '\t', stringsAsFactors = FALSE)%>%
  setNames(c("Chr","Source","Evidence", "Start", "Stop","Score", "Strand", "Frame", "Attribute")) %>%
  mutate_at(vars(Attribute), str_replace_all, pattern = ("Target Motif:"), replacement = "") %>%
  mutate_at(vars(Attribute), str_replace_all, pattern = (" .*"), replacement = "") %>%
  mutate(Strand = ifelse(.$Strand == "+", "positive", "negative")) %>%
  mutate(Attribute = paste0(.$Attribute, "_", .$Strand))

loc_df <- data.frame()
# Loop over each chromosome in a GFF3 file and identify tandem TEs
# Take the coordinate of the start and stop, save to a dataframe
for (k in unique(gff_te$Chr)){
  
  test <- subset(gff_te, Chr == k)
  print(paste0("Beginning chromosome ", k))
  
  for (l1 in 1:nrow(test)){
    
    # Modulus operation
    if(l1 %% 1000==0) {
      # Print on the screen some message
      cat(paste0("iteration: ", l1, "\n"))
    }
    
    tmp1 <- test[l1, "Attribute"]
    l2 <- l1 +1
    tmp2 <- test[l2, "Attribute"]
    
    # Check if the two LTRs are in the same family
    tryCatch({
      if (tmp1 == tmp2){
        loc_tmp <- test[l1:l2,]
        loc_tmp <- data.frame(
          loc_tmp[1, "Chr"],
          loc_tmp[1, "Start"],
          loc_tmp[2, "Stop"],
          loc_tmp[1, "Attribute"]
        ) %>%
          setNames(., nm = c("Chr", "Start", "Stop", "TE"))
        
        loc_df <- rbind(loc_df, loc_tmp)
      }
    }, error = function(e){})
  }
}

# Take the dataframe of GFF coordinates from the above command, find overlapping intervals
# Using Genomic Ranges respects both +/- strand designations and chromosome boundaries
x <- loc_df
x <- x %>%
  mutate(Strand = ifelse(grepl("_positive", .$TE), "+", "-")) %>%
  mutate(Range = paste0(.$Start, "-", .$Stop)) %>%
  mutate(TE = gsub("_positive|_negative", "", .$TE)) %>%
  #mutate(Chr = gsub("^.*_", "", .$Chr)) %>%
  dplyr::select(., c(Chr, Range, Strand, TE))

gr <- GRanges(
  seqnames = x$Chr,
  ranges = IRanges(x$Range),
  strand = x$Strand,
  te = x$TE)

x$Group <- subjectHits(GenomicRanges::findOverlaps(gr, reduce(gr)))
x <- separate(x, Range, into = c("Start", "Stop"), sep = "-")

new_df <- data.frame()
for (i in unique(x$Group)){
  
  # Modulus operation
  if(i %% 1000==0) {
    # Print on the screen some message
    cat(paste0("iteration: ", i, "\n"))
  }
  
  tmp <- x[x$Group == i,]
  tmp <- data.frame(
    Chr = tmp[1, "Chr"],
    Start = tmp[1, "Start"],
    Stop = tmp[nrow(tmp), "Stop"],
    TE = tmp[1, "TE"]
  )
  
  new_df <- rbind(new_df, tmp)
}

write.table(new_df, "~/Downloads/tandem_te_merged.bed", sep = '\t', quote = FALSE, row.names = FALSE)

