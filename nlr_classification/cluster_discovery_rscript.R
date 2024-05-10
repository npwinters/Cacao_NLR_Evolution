# Written by Noah Winters, Oct 2020
# This identified NLR clusters on each chromosome based on a window size cutoff (in genes)
# All it takes as input is:
# (1) Genotype ID (as it appears in file), 
# (2) Window size, 
# (3) List of NLR genes and their positions, tab delimited -- includes Chr, Start, Stop, Gene ID 
# (4) List of all genes and their positions, tab delimited (including NLRs)--includes Chr, Start, Stop, Gene ID 
# (5) Orthogroup file, like those from PlantTribes

list.of.packages <- c("tidyverse")
new.packages <- list.of.packages[!(list.of.packages %in% (.packages()))]
if(length(new.packages)>0){ 
  print("Loading missing packages")
  library(new.packages, character.only = TRUE, verbose = FALSE, quietly = TRUE, warn.conflicts = FALSE)
}

print ("-----------------------------------------------------------------------------------------------------------------------")

print ("To run: rscript cluster_discovery_rscript.R genotype window_size nlr_position_file gene_position_file orthogroup_file")

print ("-----------------------------------------------------------------------------------------------------------------------")

main <- function(){
  ##################################
  ###                            ###
  ###     CREATE DATAFRAMES      ###
  ###                            ###
  ##################################
  args <- commandArgs(trailingOnly = TRUE)
  genotype <- args[1]
  window_size <- args[2]
  nlr_file <- args[3]
  gene_file <- args[4]
  ortho_file <- args[5]
  
  # genotype <- "NA-807"
  # window_size <- 8
  # nlr_file <- "~/GradSchool/Research/Projects/NSF_PlantGenome/NLR_Comparative_Genomics/ideogram_files/NA-807_nlrPos.txt"
  # gene_file <- "~/GradSchool/Research/Projects/NSF_PlantGenome/NLR_Comparative_Genomics/ideogram_files/NA-807_genePos.txt"
  # ortho_file <- "~/GradSchool/Research/Projects/NSF_PlantGenome/NLR_Comparative_Genomics/evol_analyses/nlrGeneFamilyClassification_Aug2021/proteins.both.37Gv1.0.bestOrthos"

  # What genotype 
  gntype <- as.character(genotype)

  # What window size should be used, i.e. how are NLR clusters being defined 
  wind_size <- as.numeric(window_size)
  
  # NLRs
  if (gntype == "Matina"){
    # Table of identified NLRs
    nlr <- read.table(nlr_file, stringsAsFactors = FALSE, header = TRUE) %>%
      dplyr::select(Chr, NLR_Gene) %>% 
      setNames(., nm = c("Chr", "gene_id")) %>%
      filter(grepl("Chr[1-9].*", Chr) | grepl("scaffold18", Chr) | grepl("scaffold40", Chr)) %>%
      arrange(Chr, gene_id) %>%
      group_by(Chr) %>%
      group_split(.) %>% # Convert to list
      lapply(., as.data.frame)
    
    # What chromosomes are represented in the NLRs
    rep_chr <- read.table(nlr_file, stringsAsFactors = FALSE, header = TRUE) %>%
      dplyr::select(Chr, NLR_Gene) %>% 
      setNames(., nm = c("Chr", "gene_id")) %>% 
      filter(grepl("Chr[1-9].*", Chr) | grepl("scaffold18", Chr) | grepl("scaffold40", Chr)) %>%
      dplyr::select(Chr) %>% # Remove Chr0 and keep the two scaffolds with NLRs on them
      unique() %>% 
      deframe()
  } else {
    # Table of identified NLRs
    nlr <- read.table(nlr_file, stringsAsFactors = FALSE, header = TRUE) %>%
      dplyr::select(Chr, NLR_Gene) %>% 
      setNames(., nm = c("Chr", "gene_id")) %>%
      .[grep("[C|c]hr[1-9].*", .$Chr),] %>%
      arrange(Chr, gene_id) %>%
      group_by(Chr) %>%
      group_split(.) %>% # Convert to list
      lapply(., as.data.frame)
    
    # What chromosomes are represented in the NLRs
    rep_chr <- read.table(nlr_file, stringsAsFactors = FALSE, header = TRUE) %>%
      dplyr::select(Chr, NLR_Gene) %>% 
      setNames(., nm = c("Chr", "gene_id")) %>% 
      .[grep("[C|c]hr[1-9].*", .$Chr),] %>% # Remove Chr0 
      dplyr::select(Chr) %>% 
      unique() %>% 
      deframe()
  }
  
  # Genes (inclusive of NLRs)
  if (gntype == "Matina"){
    # Table of all annotated genes, inclusive of NLRs
    genes <- read.table(gene_file, stringsAsFactors = FALSE, header = TRUE) %>% 
      dplyr::select(Chr, Gene) %>% 
      setNames(., nm = c("Chr", "gene_id")) %>%
      .[.$Chr %in% rep_chr,] %>%
      arrange(Chr, gene_id) %>%
      group_by(Chr) %>%
      group_split(.) %>% # Convert to list
      lapply(., as.data.frame)
  } else {
    # Table of all annotated genes, inclusive of NLRs
    genes <- read.table(gene_file, stringsAsFactors = FALSE, header = TRUE) %>% 
      dplyr::select(Chr, Gene) %>% 
      # .[grep("[C|c]hr[1-9].*", .$Chr),] %>%
      .[.$Chr %in% rep_chr,] %>%
      setNames(., nm = c("Chr", "gene_id")) %>%
      arrange(Chr, gene_id) %>%
      group_by(Chr) %>%
      group_split(.) %>% # Convert to list
      lapply(., as.data.frame)
  }
  
  # Table of orthogroup classifications for each NLR
  orthogroups <- read.table(ortho_file, header = TRUE, stringsAsFactors = FALSE, sep = '\t') %>%
    setNames(., nm = c("gene_id", "orthogroup_id"))
  
  
  if (gntype == "Criollo"){
    orthogroups <- 
      orthogroups %>%
      .[grep("^Tc.*", .$gene_id),] %>% 
      mutate(gene_id = gsub("\\.[0-9].*", "", gene_id)) %>% 
      mutate(gene_id = gsub("_p", "_g", gene_id))
  } else if (gntype == "Matina"){
    orthogroups <- 
      orthogroups %>%
      .[grep("^Thecc.*", .$gene_id),] %>% 
      mutate(gene_id = gsub("\\.1\\.p", "", gene_id))
  } else {
    orthogroups <- 
      orthogroups %>%
      .[grep(gntype, .$gene_id),] %>% 
      mutate(gene_id = gsub("\\.[0-9].*", "", gene_id))
  }
  
  # Make sure all data frames are evenly divisible by window frame size, add col of window designations, merge with NLRs
  for (ele in 1:length(genes)){
    if (all(floor(nrow(genes[[ele]]) / wind_size) == nrow(genes[[ele]]) / wind_size, na.rm = TRUE) == FALSE){
      bot <- floor(nrow(genes[[ele]]) / wind_size)
      top <- bot + 1 
      n_rows <- top * wind_size
      add_nrows <- (top * wind_size) - nrow(genes[[ele]])
      genes[[ele]][nrow(genes[[ele]]) + add_nrows,] <- NA # Add extra NA rows 
      genes[[ele]] <- mutate(genes[[ele]], windows = rep(1:floor(nrow(genes[[ele]])/wind_size), each = wind_size))  # Add col of window IDs
      genes[[ele]] <- mutate(genes[[ele]], nlr = ifelse(genes[[ele]]$gene_id %in% nlr[[ele]]$gene_id, 1, 0)) # Which genes are NLRs?
    } else {
      #print (paste0("Is divisible by ", wind_size))
      genes[[ele]] <- mutate(genes[[ele]], windows = rep(1:floor(nrow(genes[[ele]])/wind_size), each = wind_size))
      genes[[ele]] <- mutate(genes[[ele]], nlr = ifelse(genes[[ele]]$gene_id %in% nlr[[ele]]$gene_id, 1, 0))
    }
  }
  
  genes_df <- genes
  
  ###################################
  ###                            ###
  ###    CLUSTERING ALGORITHM    ###
  ###                            ###
  ###################################
  
  # Is the data frame divisible by the window size?
  # If not, make it divisible by window size
  # if (all(floor(nrow(genes) / wind_size) == nrow(genes) / wind_size, na.rm = TRUE) == FALSE){
  #   bot <- floor(nrow(genes) / wind_size)
  #   top <- bot + 1 
  #   n_rows <- top * wind_size
  #   add_nrows <- (top * wind_size) - nrow(genes)
  #   dummy_df <- data.frame(replicate(add_nrows, c("NA", "NA") , simplify = FALSE)) %>%
  #     setNames(., nm = colnames(genes))
  #   genes <- rbind(genes,dummy_df)
  #   genes$nlr <- as.numeric(genes$nlr)
  # } else {
  #   print (paste0("Is divisible by ", wind_size))
  #   genes$nlr <- as.numeric(genes$nlr)
  # }
  
  # Create a column containing the window designation for each gene
  # genes <- genes %>%
  #   mutate(windows = rep(1:floor(nrow(genes)/wind_size), each = wind_size))
  
  # Summarize the number of NLRs in each window
  for (ele in 1:length(genes)) {
    genes[[ele]] <- 
      genes[[ele]] %>%
      group_by(windows) %>% 
      summarise(freq = sum(nlr, na.rm = TRUE)) %>%
      deframe(.)
  }
  
  # genes_freq <- genes %>%
  #   group_by(windows) %>% 
  #   summarise(freq = sum(nlr, na.rm = TRUE)) %>%
  #   deframe(.)
  
  # Merge windows when each contain at least 1 NLR
  list_clust_freq <- list()
  for (ele in 1:length(genes)) {
    genes_freq <- unlist(genes[ele])
    names(genes_freq) <- names(genes[[ele]])
    len <- length(genes_freq)
    #print (len)
    tryCatch({ # Error thrown when last window tries to merge with nothing
      for (i in seq_along(genes_freq)){
        if (genes_freq[i] > 0 && genes_freq[i+1] > 0) {
          genes_freq[i] = genes_freq[i] + genes_freq[i+1] # Merge adjacent windows if they both contain NLRs 
          names(genes_freq)[i] = paste0(names(genes_freq[i]),",", names(genes_freq[i+1]))
          genes_freq = genes_freq[-(i+1)] # Remove i + 1 element after i and i + 1 have been consolidated
        } else {
          genes_freq[i] = genes_freq[i] 
          names(genes_freq)[i] = names(genes_freq)[i]
        }
      }
    }, error = function(e){}) # Need dummy function to handle error thrown in if else statement
    list_clust_freq <- append(list_clust_freq, list(genes_freq)) # Append consolidated lists to empty list
  }
  
  
  # Assign clusters
  # clusters <- list_clust_freq
  # for (i in 1:length(clusters)) {
  #   names(clusters[[i]]) <- rep(paste0("C",1:length(clusters[[i]])), by = 1)
  # }
  
  # Convert list of named vectors to data frames and assign clusters based on merged windows 
  clusters <- list_clust_freq
  
  for (i in 1:length(clusters)){
    clusters[[i]] <- data.frame(clusters[[i]]) %>%
      setNames(., nm = "cnt") %>%
      rownames_to_column(., var = "windows") %>% 
      rownames_to_column(., var = "clusters") %>%
      mutate(., clusters = paste0("C", clusters)) %>%
      separate_rows(., windows, sep = ",", convert = TRUE) 
    clusters[[i]] <- left_join(genes_df[[i]], clusters[[i]], by = "windows") %>% 
      mutate(., windows = paste0("W", windows)) %>%
      mutate(., win_clust = paste0(Chr, "_", windows)) %>%
      mutate(., chr_clust = paste0(Chr, "_", clusters))}
  
  # Create dataframe indicating which window/cluster/ortho each gene belongs to
  clusters_df <- do.call(rbind, clusters)
  clusters_df_nlr <- clusters_df %>%
    subset(nlr > 0) %>%
    inner_join(., orthogroups, by = "gene_id") %>%
    group_by(chr_clust) %>%
    mutate(sum_nlr = sum(nlr)) %>%
    subset(sum_nlr > 1)
  
  # Write out file containing window/cluster membership for each gene
  write.csv(clusters_df_nlr, paste0(genotype,"_nlr_cluster_membership.csv"), 
            quote = FALSE, 
            row.names = FALSE)
  
  # Create ideogram file
  ortho_tally <- 
    clusters_df_nlr %>%
    group_by(chr_clust) %>%
    summarize(unique(orthogroup_id)) %>%
    tally(name = "orthogroup_sum")
  
  gene_pos <- read.table(gene_file, stringsAsFactors = FALSE, header = TRUE)  %>%
    dplyr::select(Gene, Start, Stop) %>%
    setNames(., nm = c("gene_id", "Start", "Stop"))
  
  ortho_cluster_pos <-
    clusters_df_nlr %>%
    inner_join(., gene_pos, by = "gene_id") %>%
    group_by(chr_clust) %>%
    summarize(clust_start = min(Start),
              clust_stop = max(Stop))  %>%
    inner_join(., clusters_df_nlr, by = "chr_clust") %>%
    inner_join(., ortho_tally, by = "chr_clust") %>%
    dplyr::select(Chr, clust_start, clust_stop, orthogroup_sum) %>%
    add_column(Type = ifelse(.$orthogroup_sum > 1, "Heterogeneous", "Homogeneous"), .before = 1) %>%
    add_column(Shape = ifelse(.$Type == "Heterogeneous", "circle", "triangle"), .after = 1) %>%
    dplyr::select(-orthogroup_sum) %>%
    setNames(nm = c("Type", "Shape", "Chr", "Start", "End"))
  
  write.csv(ortho_cluster_pos, paste0(genotype,"_cluster_pos_perChr_ideogram.csv"),
            quote = FALSE,
            row.names = FALSE)
  
  # clusters <- data.frame(clusters) %>%
  #   setNames(., nm = "windows") %>%
  #   rownames_to_column(., var = "clusters") %>%
  #   separate_rows(., windows, sep = ",", convert = TRUE) %>%
  #   left_join(genes,., by = "windows")
  
  # How many NLRs are in each cluster?
  # Cluster here defined as 2 or more genes that are within 9 genes of one another 
  nlr_clust_sum <- clusters 
  
  for (i in 1:length(nlr_clust_sum)){
    nlr_clust_sum[[i]] <- nlr_clust_sum[[i]] %>%
      group_by(clusters, Chr) %>%
      summarise(sum_nlr = sum(nlr)) %>%
      mutate(., chr_clust = paste0(Chr, "_", clusters)) %>%
      subset(sum_nlr > 1)}
  
  # Dataframe containing all clusters with > 1 NLR 
  nlr_clust_sum_df <- do.call(rbind, nlr_clust_sum)
  
  # Write out file containing NLR sums per cluster
  write.csv(nlr_clust_sum_df, paste0(genotype,"_nlr_sum_clusters_perChr.csv"), 
            quote = FALSE, 
            row.names = FALSE)
  
  # On which chromosomes do we find most of the NLRS?
  cluster_totals <- nlr_clust_sum_df %>% 
    group_by(Chr) %>% 
    summarise(length(clusters)) %>% 
    setNames(., nm = c("Chr", "No. NLR Clusters"))
  
  # Write out file containing # NLR clusters per chromosome
  write.csv(cluster_totals, paste0(genotype,"_cluster_totals_perChr.csv"), 
            quote = FALSE, 
            row.names = FALSE)
  
  print ("---------------------------------------------------------------------------------")
  print (paste0(gntype," NLR Clusters:"))
  print (cluster_totals)
  print ("---------------------------------------------------------------------------------")
}

main()



