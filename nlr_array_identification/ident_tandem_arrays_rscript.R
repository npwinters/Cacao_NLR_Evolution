#!/usr/bin/env Rscript
# Written by Noah Winters, Sept 2021
# This script identifies tandem arrays on each chromosome based on output from MCScanX

# Load packages 
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("tidyverse"))

##################################
###                            ###
###      Set up argparse       ###
###                            ###
##################################
# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-g", "--genotype", type="character", default="_outFile", 
                    help="Genotype of interest [default %(default)s]")
parser$add_argument("tandem_file", nargs = 1, 
                    help = "Tab-delimited tandem array file from parse_tandem_arrays.sh. [default %(default)s]")
parser$add_argument("gene_file", nargs = '?', 
                    help="File containing gene IDs of interest (such as those from from classify_nlr.sh). [default %(default)s]")

# Get command line options, if help option encountered print help and exit,
# Otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()

###############################
###                         ###
###     PRE-PROCESSING      ###
###                         ###
###############################
# Arguments
gntype <- as.character(args$genotype)

# Read in file of parsed tandem arrays and NLR 
genes_df <- read.table(args$tandem_file, header = TRUE) %>% 
  mutate(Dup_Type = ifelse(Dup_Type == "Tandem", "Tandem", .$Position)) # Need the Position col to be == row position

if (length(args$gene_file) > 0){
  args$gene_file
  nlr <- read.table(args$gene_file) %>%
    setNames(nm = "Gene") %>%
    mutate(Gene = gsub("\\.1\\.p", "", Gene)) %>% 
    mutate(Gene = gsub("_p", "_g", Gene)) %>%
    mutate(Gene = gsub("\\.1", "", Gene))
}

# Which chromosome are we working with?
if (gntype == "[C|c]riollo"){
  Chr <- as.character(genes_df[1,"Gene"]) %>%
    gsub("v2_.*", "", .) %>%
    gsub("Tc", "Criollo_Chr", .)
} else if (gntype == "[M|m]atina"){
  Chr <- as.character(genes_df[1,"Gene"]) %>%
    gsub("G.*", "", .) %>%
    gsub("Thecc\\.", "Matina_Chr")
} else if (gntype != "[M|m]atina" | gntype != "[C|c]riollo"){
  Chr <- as.character(genes_df[1,"Gene"]) %>%
    gsub("v1_.*", "", .)
} 

##################################
###                            ###
###     CLUSTERING ARRAYS      ###
###                            ###
##################################
# Discover tandem arrays using RLE to find repeated cases of "tandem" 
# Note: The tandem arrays identified in MCScanX's HTML files differ slightly from the gene_dup file
# This seems to mainly be a problem for cases when a tandem dup is labeled segmental/WGD

# Identifies runs of tandem duplicates that are > 3, find their index values
runs = rle(as.character(genes_df$Dup_Type))
tandems <- which(runs$values == "Tandem" & runs$lengths >= 3)

# If chromosome _begins_ with a tandem array, adjust start position to first gene
if ( tandems[1] == 1 ){
  start <- as.numeric((runs$values[tandems-1])) + 1 
  start <- append(start, 1, after = 0)
} else {
  start <- as.numeric((runs$values[tandems-1])) + 1 
}

# If chromosome _ends_ with a tandem array, adjust end position to last gene
if ( any(is.na(as.numeric((runs$values[tandems+1])))) ){
  end <- as.numeric((runs$values[tandems+1])) - 1
  end[length(end)] <- nrow(genes_df)
} else {
  end <- as.numeric((runs$values[tandems+1])) - 1
}

# Subset the dataframe by the runs of tandem duplicates
arrays <- list()
if (length(start) > 0){
  for ( i in 1:length(start)){
    array_number <- paste0(Chr, "_array", i)
    tmp <- as.character(genes_df[start[i]:end[i],]$Gene)
    arrays[[array_number]] <- tmp
    } 
  } else {
  cat("There are no tandem arrays on this chromosome.\n")
} 

# Create dataframe
if (length(arrays) > 0){
  if (length(args$gene_file > 0)){
    array_df <- plyr::ldply(arrays, rbind) %>%
      rename("array_id" = ".id") %>%
      group_by(array_id) %>%
      pivot_longer(cols = -array_id, values_to = "Gene") %>% 
      na.omit() %>%
      dplyr::select(-name) %>% 
      add_column(is_nlr = ifelse(.$Gene %in% nlr$Gene, 
                                 yes = "1", 
                                 no = "0"))
  } else {
    array_df <- plyr::ldply(arrays, rbind) %>%
      rename("array_id" = ".id") %>%
      group_by(array_id) %>%
      pivot_longer(cols = -array_id, values_to = "Gene") %>% 
      na.omit() %>%
      dplyr::select(-name)
  }
  
  # Write output files
  outfile <- paste0(Chr, "_tandemArrays.txt")
  write.table(array_df, 
              file = outfile, 
              sep = '\t', 
              quote = FALSE, 
              row.names = FALSE)
  
}

  
  
  

  










  
