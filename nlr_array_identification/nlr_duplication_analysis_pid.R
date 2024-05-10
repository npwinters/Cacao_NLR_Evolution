library(tidyverse)
library(patchwork)
library(gridExtra)
library(broom)
library(ggpmisc)
library(equatiomatic)

setwd("/Users/noah/GradSchool/Research/Projects/NSF_PlantGenome/NLR_Comparative_Genomics/")
pid_file <- "evol_analyses/phylogenomics/all_10x_pep_10aug2021_nbarc_rpw8_tir_hmm_nlrparser.pid"
pid <- read.table(pid_file, header = FALSE, row.names = 1, skip = 1)
colnames(pid) <- rownames(pid)


gntype <- c("Criollo", "Matina", "SCA-6", "SPEC-54-1", "CCN-51",
            "GU-257E", "ICS-1", "IMC-105",  "NA-246", "NA-807",
            "Pound-7")
mod_sums_total <- data.frame(
  Chr = character(),
  Rsquared = numeric(), 
  Pvalue = numeric()
)
grid_plot_list <- list()

for (i in gntype){
  
  # Read in percent identity matrix
  gene_position_file = paste0("ideogram_files/", i, "_genePos.txt")
  if (i == "Criollo"){
    gene_pos <- read.table(gene_position_file, header = TRUE) %>%
      mutate(Chr = gsub("chr", "Criollo_Chr", Chr))
  } else {
    gene_pos <- read.table(gene_position_file, header = TRUE)
  }
  
  # Subset the percent identity matrix by genotype
  if (i == "Matina"){
    pid_gntype <- pid[grep("Thecc", rownames(pid)),] %>% 
      .[,grep("Thecc", colnames(.))]
  } else if(i == "Criollo"){
    pid_gntype <- pid[grep("Tc", rownames(pid)),] %>% 
      .[,grep("Tc", colnames(.))]
  } else {
    pid_gntype <- pid[grep(i, rownames(pid)),] %>% 
      .[,grep(i, colnames(.))]
  }
  
  
  # Calculate distance between genes 
  tPID <- t(combn(colnames(pid_gntype), 2))
  tPID <- data.frame(tPID, dist=pid[tPID]) %>% 
    setNames(., nm = c("Gene", "Gene2", "PID")) %>%
    mutate(Gene = gsub("\\.1\\.p", "", .$Gene)) %>%
    mutate(Gene2 = gsub("\\.1\\.p", "", Gene2)) %>%
    mutate(Gene = gsub("\\.1", "", .$Gene)) %>%
    mutate(Gene2 = gsub("\\.1", "", Gene2)) %>%
    mutate(Gene = gsub("_p", "_g", .$Gene)) %>%
    mutate(Gene2 = gsub("_p", "_g", Gene2)) %>%
    inner_join(gene_pos, by = "Gene") %>%
    rename("Gene1_Start"=Start,
           "Gene1_Stop"=Stop,
           "Gene1" = Gene, 
           "Gene" = Gene2) %>%
    inner_join(rename(gene_pos, "Gene2_Chr" = Chr), by = "Gene") %>%
    rename("Gene2_Start"=Start,
           "Gene2_Stop"=Stop,
           "Gene2" = Gene, 
           "Gene1_Chr" = Chr) %>% 
    group_by(Gene1) %>%
    mutate(maxG1 = max(c(Gene1_Start, Gene1_Stop))) %>%
    group_by(Gene2) %>%
    mutate(minG2 = min(c(Gene2_Start, Gene2_Stop))) %>% 
    ungroup() %>% 
    mutate(dist = minG2 - maxG1)
  
  #  Make sure we are comparing across a single chromosome
  # Remove negative values, i.e. genes are are overlapping one another on opposite strands
  tPID <- tPID[tPID$Gene1_Chr == tPID$Gene2_Chr,] %>%
      filter(dist > 0)
  
  # Initiate empty data frame
  mod_sums <- data.frame(
    Chr = character(),
    Rsquared = numeric(), 
    Pvalue = numeric()
  )
  
  # Calculate correlation between distance and percent identity
  for (k in as.character(unique(tPID$Gene1_Chr))){
    tmp <- subset(tPID, Gene1_Chr == k)
    if (nrow(tmp) > 1){
      corr <- cor(tmp$PID, log2(tmp$dist), method = "spearman")
      summod <- summary(lm(tmp$PID ~ log2(tmp$dist)))
      summod <- data.frame(cbind(Chr = k, 
                                 Spearman.Corr = round(corr, 5),
                                 R.squared = round(summod$r.squared, 5), 
                                 P.value = round(summod$coefficients[2,4], 5)))
      mod_sums <- rbind(mod_sums, summod)
    } else{
      summod <- data.frame(cbind(Chr = k, 
                                 Spearman.Corr = NA,
                                 R.squared = NA, 
                                 P.value = NA))
      mod_sums <- rbind(mod_sums, summod)
      
    }
  }
  
  mod_sums_total <- rbind(mod_sums_total, mod_sums) 
  
  p1 <-
    tPID %>%
    mutate(Gene1_Chr = factor(Gene1_Chr, 
                              levels = c(paste0(i, "_Chr1"), paste0(i, "_Chr2"), paste0(i, "_Chr3"), paste0(i, "_Chr4"), paste0(i, "_Chr5"),
                                         paste0(i, "_Chr6"), paste0(i, "_Chr7"), paste0(i, "_Chr8"), paste0(i, "_Chr9"), paste0(i, "_Chr10")),
                              labels = c(paste0(i, "_Chr1"), paste0(i, "_Chr2"), paste0(i, "_Chr3"), paste0(i, "_Chr4"), paste0(i, "_Chr5"),
                                         paste0(i, "_Chr6"), paste0(i, "_Chr7"), paste0(i, "_Chr8"), paste0(i, "_Chr9"), paste0(i, "_Chr10")))) %>% 
    ggplot(aes(x = PID, y = log2(dist))) + 
    geom_point(alpha = 0.25, color = "black") + 
    geom_smooth(method = "lm", color = "red", linetype = "dotted", size = 1, fullrange = FALSE) +
    facet_wrap(~Gene1_Chr) + 
    theme_minimal_grid() + 
    labs(x = "Percent Identity", y = "log2[distance in bp]") +
    theme(
      strip.background = element_rect(color = "lightgrey")
    )
  
  
  mod_sums <- 
      mod_sums %>%
      mutate(Chr.num = str_replace_all(as.character(.$Chr), c("chrUn_random" = "11", "^[A-Z].*_Chr" = "", "^[A-Z].*_scaffold"=""))) %>%
      mutate(Chr.num = as.numeric(Chr.num)) %>%
      arrange(Chr.num) %>%
      dplyr::select(c(Chr, Spearman.Corr, R.squared, P.value))
  
  
  p1 <- ggplotGrob(p1)
  p2 <- tableGrob(mod_sums, rows = NULL)
  
  
  grid_plot <- grid.arrange(p1, p2)
  
  grid_plot_list[[i]] <- grid_plot
  
}

pdf("figures/NLR_distancePIDCorr_allGenotypes.pdf", height = 9)
lapply(grid_plot_list, grid.arrange)
dev.off()


##################
#### Analyses ####
##################

# Does not look like Distance x PID correlates with and dispersal type
mod_sums_nlr <- 
  nlr %>% 
  add_column(Genotype_Chr = paste0(.$Genotype, "_Chr", .$Chr)) %>% 
  group_by(Genotype_Chr, Dup_Type) %>%
  tally(name = "Count") %>%
  rename("Chr" = Genotype_Chr) %>%
  inner_join(mod_sums_total, by = "Chr")

mod_sums_nlr <- na.omit(mod_sums_nlr) %>%
  pivot_wider(names_from = Dup_Type, values_from = Count) %>% 
  mutate(Spearman.Corr = as.numeric(as.character(Spearman.Corr)))

mod_sums_nlr[is.na(mod_sums_nlr)] <- 0
mod_sums_nlr$total_nlr <- rowSums(mod_sums_nlr[,c("tandem", "dispersed", "proximal", "segmental", "singleton")])

pid_mod <- aov(Spearman.Corr ~ dispersed + tandem + proximal + segmental + singleton, mod_sums_nlr)
tidy(pid_mod)

plot(mod_sums_nlr$Spearman.Corr, mod_sums_nlr$tandem)


# Do dispersed duplicates drive heterogeneous cluster formation?
cluster_pos <- read.csv("evol_analyses/nlr_clusters/allGntypes_cluster_pos_perChr_ideogram.csv")

cluster_pos <- cluster_pos %>% 
  group_by(Chr, Type) %>% 
  tally() %>% 
  mutate(combChr = gsub("^[A-Z].*_Chr", "", Chr)) %>%
  inner_join(mod_sums_nlr, by = "Chr") %>%
  pivot_wider(names_from = Type, values_from = n) 
cluster_pos[is.na(cluster_pos)] <- 0
cluster_pos <- 
  cluster_pos %>%
  mutate(propHet = Heterogeneous / (Heterogeneous + Homogeneous),
         propDisp = dispersed / (singleton + dispersed + tandem + proximal + segmental))

# Make the plot
p1 <-
  cluster_pos %>%
  rename("total" = total_nlr) %>%
  pivot_longer(cols = c(segmental, singleton, proximal, dispersed, tandem, total),
               names_to = "Dup_Type", values_to="Count") %>% 
  mutate(Dup_Type = factor(Dup_Type, 
                           levels = c("singleton", "dispersed", "proximal", "tandem", "segmental", "total"),
                           labels = c("singleton", "dispersed", "proximal", "tandem", "segmental", "total"))) %>%
  subset(Dup_Type != "singleton") %>%
  ggplot(aes(x = Count, y = Heterogeneous)) + 
  geom_point(show.legend = FALSE) + 
  geom_smooth(formula = y ~ x, method = "lm", linetype = "dotted") +
  labs(x = "# NLR Genes", y = "# Heterogeneous Clusters") +
  theme_minimal_grid() +
  facet_wrap(~Dup_Type, ncol = 6, nrow = 1, scales = "free") + 
  theme(
    strip.background = element_rect(color = "darkgrey")
  )

# Construct the model
mod1 <- summary(lm(Heterogeneous ~  proximal + tandem + segmental + dispersed, cluster_pos)) %>%
  .$coefficients %>%
  data.frame() %>%
  setNames(nm = c("Estimate", "Std Error", "T.value", "P.value")) %>%
  round(., digits = 5)

p1 <- ggplotGrob(p1)
p2 <- tableGrob(mod1)

pdf("figures/NLR_dupType_hetCluster_corr.pdf", height = 6, width = 9)
grid.arrange(p1, p2, padding =2)
dev.off()

# Do dispersed duplicates drive heterogeneous cluster formation--unique orthogroups?
allClust_mem <- read.csv("evol_analyses/nlr_clusters/allGntypes_nlr_cluster_membership.csv", header = TRUE)


allClust_mem <- 
  allClust_mem %>%
  dplyr::select(c(Chr, Window_Start, Window_Stop, Orthogroup)) %>%
  group_by(Chr) %>% 
  summarise(unique_orthos = n_distinct(Orthogroup)) %>%
  inner_join(cluster_pos, by = "Chr")

#pdf("figures/NLR_dupType_uniqOrtho_corr.pdf", height = 2.5, width = 9)
p1 <-
  allClust_mem %>%
  rename("total" = total_nlr) %>%
  pivot_longer(cols = c(segmental, singleton, proximal, dispersed, tandem, total),
               names_to = "Dup_Type", values_to="Gene_Count") %>% 
  mutate(Dup_Type = factor(Dup_Type, 
                           levels = c("singleton", "dispersed", "proximal", "tandem", "segmental", "total"),
                           labels = c("singleton", "dispersed", "proximal", "tandem", "segmental", "total"))) %>%
  subset(Dup_Type != "singleton") %>%
  ggplot(aes(x = Gene_Count, y = unique_orthos)) + 
  geom_point(show.legend = FALSE) + 
  geom_smooth(formula = y ~ x, method = "lm", linetype = "dotted") +
  labs(x = "# NLR Genes", y = "# Unique Orthogroups") +
  theme_minimal_grid() +
  facet_wrap(~Dup_Type, ncol = 5, nrow = 1, scales = "free") + 
  theme(
    strip.background = element_rect(color = "darkgrey")
  ) 
#dev.off()


# Construct the model
mod2 <- summary(lm(unique_orthos ~  proximal + tandem + segmental + dispersed, allClust_mem)) %>%
  .$coefficients %>%
  data.frame() %>%
  setNames(nm = c("Estimate", "Std Error", "T.value", "P.value")) %>% 
  round(., digits = 5)

p1 <- ggplotGrob(p1)
p2 <- tableGrob(mod2)

pdf("figures/NLR_dupType_uniqOrtho_corr.pdf", height = 6, width = 9)
grid.arrange(p1, eq, p2, padding = 2)
dev.off()











