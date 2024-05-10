# List of packages
packages <- c("tidyverse", "cowplot")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
lapply(packages, function(i){
  suppressPackageStartupMessages(library(i, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE))
})

# suppressPackageStartupMessages(
#   suppressWarnings(
#     lapply(list("argparse", "tidyverse", "cowplot"), library)
#     ))

print ("------------------------------------------------------------------------------------")

print ("This script categorizes NLRs into groups according to their domain structure")

print ("Usage: bash categorize_nlr.sh <arguments>")

print ("Arguments: PFAM_OUTPUT.TSV NLRPARSER_CC_ONLY.txt OUTPUT_PREFIX")

print ("------------------------------------------------------------------------------------")

main <- function(){
  
  #################
  ### Arguments ###
  #################
  args <- commandArgs(trailingOnly = TRUE)
  pfam_file <- args[1]
  parser_cc <- args[2]
  output_prefix <- args[3]
  output_prefix <- as.character(output_prefix)
  
  ####################
  ### Read in data ###
  ####################
  
  # Read in nlr-parser cc only table
  parser_cnl <- read.table(parser_cc, header = TRUE)
  
  # Read in and parse domain table
  domTable <- 
    read.table(pfam_file, sep = '\t', header = FALSE, fill = TRUE, stringsAsFactors = FALSE) %>%
    setNames(., nm = c("Gene", "Protein_accession",
                       "Seq_length", "Analysis", 
                       "Sig_accession", "Sig_description", 
                       "Start_location", "Stop_location", 
                       "E-Value", "Status", 
                       "Date", "InterPro_accession", 
                       "InterPro_description")) %>%
    add_column(., Genotype = gsub("_Chr.*", "", .$Gene) , .after = 1) %>%
    mutate(Genotype = ifelse(grepl("Tc.*", Genotype), "Criollo", Genotype)) %>%
    mutate(Genotype = ifelse(grepl("Thecc.*", Genotype), "Matina", Genotype)) %>%
    add_column(Version = ifelse(grepl("Criollo", .$Genotype), "v2.0",
                                ifelse(grepl("Matina", .$Genotype), "v2.1", "v1.0")), .after = 2) %>% 
    full_join(parser_cnl, domTable, by = "Gene")
  
  ########################
  ### Domain diversity ###
  ########################
  
  # Calculate NLR integrated domain diversity
  div_table <- 
    domTable %>%
    subset(., Analysis == "Pfam") %>%
    mutate(InterPro_description = ifelse(InterPro_description == '', Sig_description, InterPro_description)) %>% 
    mutate(InterPro_description = ifelse(Sig_accession == 'Coil', Sig_accession, InterPro_description)) %>%
    mutate(InterPro_accession = ifelse(InterPro_accession == '', Sig_accession, InterPro_accession)) %>%
    mutate(InterPro_description = tolower(InterPro_description)) %>%
    mutate(InterPro_description = gsub("leucine[ |-]rich.*", "leucine-rich repeat", InterPro_description)) %>% 
    group_by(InterPro_description) %>%
    tally()
  
  # Remove the main NLR domains
  main_domains <- c("nb-arc","leucine-rich repeat", "rx n-terminal domain", 
                    "toll/interleukin-1 receptor homology (tir) domain", 
                    "powdery mildew resistance protein, rpw8 domain", "coil")
  
  # Subset the table 
  `%notin%` <- Negate(`%in%`)
  div_table <- div_table[div_table$InterPro_description %notin% main_domains,]
  
  # Write out the table
  write.table(div_table,
              paste0(output_prefix, "_idDiversity.txt"),
              sep = '\t',
              quote = FALSE,
              row.names = FALSE)
  
  
  ##########################
  ### NLR categorization ###
  ##########################
  
  # Categorize NLRs into TNL, CNL, NL, RNL
  # All: CC-containing = CNL, TIR-containing = TNL, RPW8-containing = RNL
  # All: Only NB-ARC, and NB-ARC + LRR = NL
  # All: Only LRR = LRR
  # Can only have a CC if you have an Rx N-term domain (PF18052), or both a Coils + NLR-Parser CC annot
  nlr_cat <- 
    domTable %>%
    mutate(InterPro_description = ifelse(InterPro_description == '', Sig_description, InterPro_description)) %>% 
    mutate(InterPro_description = ifelse(Sig_accession == 'Coil', Sig_accession, InterPro_description)) %>%
    mutate(InterPro_accession = ifelse(InterPro_accession == '', Sig_accession, InterPro_accession)) %>%
    mutate(InterPro_description = tolower(InterPro_description)) %>%
    mutate(InterPro_description = gsub("leucine[ |-]rich.*", "leucine-rich repeat", InterPro_description)) %>% 
    mutate(nlr_domains = case_when(
      str_detect(InterPro_description, regex("leucine-rich", ignore_case=TRUE)) ~ "LRR",
      str_detect(InterPro_description, regex("nb-arc", ignore_case=TRUE)) ~ "NBARC",
      str_detect(InterPro_description, regex("toll.*", ignore_case=TRUE)) ~ "TIR",
      str_detect(InterPro_description, regex("rpw8", ignore_case=TRUE)) ~ "RPW8",
      #str_detect(InterPro_description, regex("rx n-terminal", ignore_case=TRUE)) ~ "CC",
      TRUE ~ "Other"
    )) %>%
    mutate(nlr_domains = ifelse(grepl("Other", nlr_domains), InterPro_accession, nlr_domains)) %>%
    mutate(nlr_domains = ifelse(InterPro_description == "rx n-terminal domain" | Sig_accession == "Coil" & nlr_parser_class == "CNL", "CC", nlr_domains)) %>%
    arrange(Gene, Start_location)
  
  # NLR categorization + integrated domain IDs
  # nlr_cat_id <- 
  #   nlr_cat %>%
  #   subset(., Analysis == "Pfam") %>%
  #   group_by(Gene) %>%
  #   mutate(nlr_domains = gsub("IPR.*", "X", nlr_domains)) %>%
  #   dplyr::summarize(dom_abbreviation=paste(nlr_domains, collapse='-')) 
  
  
  # Calculate number of NLRs in each class, per accession
  nlr_cat <-
    nlr_cat %>%
    filter(!str_detect(nlr_domains, "IPR")) %>%
    group_by(Gene) %>%
    dplyr::summarize(dom_abbreviation=paste(nlr_domains, collapse='-')) %>%
    add_column(nlr_abbreviation = case_when(
      str_detect(.$dom_abbreviation, regex("CC", ignore_case = TRUE)) ~ "CNL",
      str_detect(.$dom_abbreviation, regex("TIR", ignore_case = TRUE)) ~ "TNL",
      str_detect(.$dom_abbreviation, regex("RPW8", ignore_case = TRUE)) ~ "RNL",
      str_detect(.$dom_abbreviation, regex("CC|TIR|RPW8|NBARC", ignore_case = TRUE), negate = TRUE) ~ "LRR",
      TRUE ~ "NL"
    )) %>% 
    add_column(., Genotype = gsub("_Chr.*", "", .$Gene) , .after = 1) %>%
    mutate(Genotype = ifelse(grepl("Tc.*", Genotype), "Criollo", Genotype)) %>%
    mutate(Genotype = ifelse(grepl("Thecc.*", Genotype), "Matina", Genotype)) %>%
    add_column(Version = ifelse(grepl("Criollo", .$Genotype), "v2.0",
                                ifelse(grepl("Matina", .$Genotype), "v2.1", "v1.0")), .after = 2) %>% 
    setNames(nm = c("Gene", "Genotype", "Version", "Architecture", "Class"))
  
  
  # Write out the table
  write.table(nlr_cat,
              paste0(output_prefix, "_nlrClass.txt"),
              sep = '\t',
              quote = FALSE,
              row.names = FALSE)
  
  # Create confidence intervals
  # sum_stats <- 
  #   nlr_cat %>%
  #   group_by(Genotype, Class) %>%
  #   tally() %>%
  #   subset(Class != "LRR") %>%
  #   summarySE(measurevar = "n", groupvars = "Class") %>%
  #   pivot_longer(cols = contains("Class"), values_to = "nlr_type")
  
  sum_stats <- 
    nlr_cat %>%
    group_by(Genotype, Class) %>%
    tally()  %>% 
    pivot_wider(names_from = Class, values_from = n) %>% 
    dplyr::select(-LRR) %>%
    add_column(Total = rowSums(.[,-1], na.rm = TRUE)) %>%
    pivot_longer(cols = c(NL, CNL, TNL, RNL, Total), names_to = "Class", values_to = "n")
  
  sum_stats_CI <-
    lapply(unique(sum_stats$Class),
           function(x){
             confint(lm(n ~ 1, sum_stats[sum_stats$Class==x,]))
           })
  
  sum_stats_CI <- 
    sum_stats_CI %>% 
    unlist(.) %>%
    matrix(.) %>%
    data.frame(.) %>% 
    setNames(., nm = "ConfInt") %>%
    add_column(Class = rep(unique(sum_stats$Class), each = 2), .before = 1) %>%
    add_column(minMax = rep(c("ymin", "ymax"), times = 5), .after = 1) %>%
    pivot_wider(names_from = minMax, values_from = ConfInt) 
  
  lbl <- 
    nlr_cat %>%
    group_by(Genotype, Class) %>%
    tally() %>%
    pivot_wider(names_from = Class, values_from = n) %>% 
    dplyr::select(-LRR) %>%
    add_column(Total = rowSums(.[,-1], na.rm = TRUE)) %>%
    pivot_longer(cols = c(NL, CNL, TNL, RNL, Total), names_to = "Class", values_to = "n") %>%
    mutate(Class = factor(Class,
                          levels = c("NL", "CNL", "TNL", "RNL", "Total"),
                          labels = c("NL", "CNL", "TNL", "RNL", "Total")))
  
  # Plot NLRs by accession and class--mean +/- 95% C.I. 
  nlr_cat %>%
    group_by(Genotype, Class) %>%
    tally() %>%
    pivot_wider(names_from = Class, values_from = n) %>% 
    dplyr::select(-LRR) %>%
    add_column(Total = rowSums(.[,-1], na.rm = TRUE)) %>%
    pivot_longer(cols = c(NL, CNL, TNL, RNL, Total), names_to = "Class", values_to = "n") %>%
    mutate(Class = factor(Class,
                          levels = c("NL", "CNL", "TNL", "RNL", "Total"),
                          labels = c("NL", "CNL", "TNL", "RNL", "Total"))) %>%
    ggplot() + 
    geom_jitter(aes(x = Class, y = n, color = Class), 
                size = 5, show.legend = TRUE, width = 0.1) + 
    stat_summary(aes(x = Class, y = n, color = Class), color = "black", shape = 5, size = 2, fun = mean) +
    labs(y = "NLR Genes per Accession", x = NULL, color = NULL) + 
    geom_errorbar(data = sum_stats_CI, aes(x = Class, ymin = ymin, ymax = ymax ), width = 0.0) + 
    theme_minimal_grid() +
    scale_color_manual(values = c("#046C9A", "#D69C4E", "#ABDDDE", "#000000", "darkgrey")) +
    theme(
      legend.position = "none"
    ) +
    geom_label_repel(data = subset(lbl, Class == "Total"),
                    aes(x = Class, y = n, label = Genotype),
                    fill = "white", min.segment.length = Inf)
  
  ggsave("NLR_Genes_perAccession_byClass.pdf", device = "pdf", width = 7, height = 7)
}
main()






            

