library(tidyverse)
library(cowplot)
library(scales)
library(performance)
library(multcompView)
library(Hmisc)
library(performance)
library(MASS)
library(emmeans)
setwd("/Users/noah/GradSchool/Research/Projects/NSF_PlantGenome/NLR_Comparative_Genomics")

#################
### Arguments ### 
#################
dup_file <- "evol_analyses/MCScanX/all_10x_pep_10aug2021_geneDupType.txt"
class_file <- "annotation/pfam_output/all_10x_pep_10aug2021_nbarc_rpw8_tir_hmm_nlrparser_nlrClass.txt"
ortho_file <- "evol_analyses/nlrGeneFamilyClassification_Aug2021/proteins.both.37Gv1.0.bestOrthos"
gene_content <- "genome_files/proteins/gene_content_perGntype.txt"
aed_scores <- "genome_files/gff/meta_assembly_annotations/metaAssembly_annot_aedScores.txt"
gene_pos_file <- "ideogram_files/allGntype_genePos.txt"

####################
### Read in data ### 
####################

# Categorization of each genes' duplication status using MCScanX (Wang et al., 2012)
# 0 = singleton
# 1 = dispersed duplicate
# 2 = proximal duplicate
# 3 = tandem duplicate
# 4 = segmental duplicate
dup <- read.table(dup_file, header = FALSE) %>%
  setNames(nm = c("Gene", "Dup_Type")) %>%
  mutate(Dup_Type = ifelse(Dup_Type == 0, "singleton", 
                           ifelse(Dup_Type == 1, "dispersed",
                                  ifelse(Dup_Type == 2, "proximal",
                                         ifelse(Dup_Type == 3, "tandem", "segmental")))))

# NLR class categorization table--genes without Pfam domains are excluded 
class <- read.table(class_file, header = TRUE)

# The orthogroup classification contains the complete set of genes predicted to be NLRs
# The NLR class categorization is a smaller subset and excludes genes without any Pfam domains 
# The difference is small, only ~30 genes
ortho <- read.table(ortho_file, sep = '\t', header = TRUE) %>%
  setNames(nm = c("Gene", "Orthogroup")) %>%
  inner_join(class, by = "Gene") %>%
  mutate(Gene = ifelse(Genotype == "Matina", gsub("\\.1\\.p", "", Gene), gsub("\\.1", "", Gene))) %>%
  mutate(Gene = gsub("_p", "_g", Gene))

# Merge data frames
nlr <- inner_join(dup, ortho, by = "Gene") %>%
  add_column(is_nlr = rep("NLR", nrow(.)))

nlr$Chr <- str_replace_all(as.character(nlr$Gene), c("Thecc\\." = "", "G[0-9].*" = "")) %>%
  str_replace_all(c("Tc" = "", "v2_[a-z].*" = "")) %>% 
  str_replace_all(c("Chr0v1" = "0", "^[A-Z].*_Chr" = "", "v1_.*" = "")) %>% 
  str_replace_all(c("^0" = "", "Un" = "0", "[A-Z].*" = "0"))

non_nlr <- anti_join(dup, ortho, by = "Gene") %>%
  add_column(Genotype = ifelse(grepl("Thecc", .$Gene), "Matina",
                               ifelse(grepl("Tc.*", .$Gene), "Criollo", gsub("_Chr.*", "", .$Gene))),
             is_nlr = rep("Non-NLR", nrow(.)))

# Function to calculate sem
sem <- function(x) sqrt(var(x)/length(x))


#########################
### NLR Summary Stats ### 
#########################
# Summarise High vs Low CNV mean and SEM
nlr %>% 
  group_by(Genotype) %>% 
  tally() %>% 
  add_column(CNV = ifelse(.$Genotype == "CCN-51" | .$Genotype == "ICS-1" |
                            .$Genotype == "GU-257E" | .$Genotype == "NA-246", "Low", "High")) %>% 
  group_by(CNV) %>% 
  summarise(mean(n),
            sem(n))

## Are high vs low significant different from one another?
nlr %>% 
  group_by(Genotype) %>% 
  tally() %>% 
  add_column(CNV = ifelse(.$Genotype == "CCN-51" | .$Genotype == "ICS-1" |
                            .$Genotype == "GU-257E" | .$Genotype == "NA-246", "Low", "High")) %>% 
  wilcox.test(n ~ CNV, .)
  

# Summarise mean and SEM for each NLR class
nlr %>% 
  group_by(Genotype, Class) %>% 
  tally() %>% 
  group_by(Class) %>% 
  summarise(mean(n),
            sem(n))
  

## Are the differences between high vs low for each class significant?
mod <- 
  nlr %>% 
  group_by(Genotype, Class) %>% 
  tally(name = "num_nlr") %>%
  add_column(CNV = ifelse(.$Genotype == "CCN-51" | .$Genotype == "ICS-1" |
                            .$Genotype == "GU-257E" | .$Genotype == "NA-246", "Low", "High"))


## Fit a GLM with a poisson distribution--good place to start for count data
glm_mod <- glm(num_nlr ~ CNV:Class, family = poisson, data = mod)

## Does a poisson dist fit the data? If p > 0.05 a poisson dist is appropriate
1 - pchisq(summary(glm_mod)$deviance, summary(glm_mod)$df.residual)

## Fit a GLM with a negative binomial distribution--poisson didn't seem to fit
glm_mod <- glm.nb(num_nlr ~ CNV:Class,  data = mod)

## Does a negative binomial dist fit the data? If p > 0.05 a negative binomial dist is appropriate
1 - pchisq(summary(glm_mod)$deviance, summary(glm_mod)$df.residual)

## Use emmeans to do posthoc comparisons
pairs(emmeans(glm_mod, ~ CNV:Class))




############
### Plot ### 
############
# NLR Duplications Relative to the Whole Genome
sum_stats <- 
  bind_rows(nlr, non_nlr) %>%
  group_by(is_nlr, Genotype, Dup_Type) %>%
  tally(name = "Count") %>%
  group_by(is_nlr, Genotype) %>%
  mutate(Proportion = Count/sum(Count)) %>%
  mutate(Dup_Type = factor(Dup_Type, levels = c("singleton","dispersed", "proximal", "tandem", "segmental"),
                           labels = c("singleton","dispersed", "proximal", "tandem", "segmental"))) %>% 
  add_column(CNV = ifelse(.$Genotype == "CCN-51" | .$Genotype == "ICS-1" |
                            .$Genotype == "GU-257E" | .$Genotype == "NA-246", "Low", "High"))

# # Code to manually add confidence intervals 
# for (i in unique(sum_stats$CNV)){ 
#   if (i == "High"){
#     tmp <- subset(sum_stats, CNV == i & is_nlr == "NLR")
#     high_ci <-
#       sapply(unique(tmp$Dup_Type),
#              function(k){
#                confint(lm(Count ~ 1, subset(tmp, Dup_Type == k)))
#              }, simplify = TRUE) %>% 
#       data.frame() %>%
#       setNames(nm = unique(tmp$Dup_Type)) %>% 
#       add_column(minMax = c("ymin", "ymax")) %>%
#       mutate(CNV = rep("High", nrow(.))) %>% 
#       pivot_longer(!c(minMax, CNV), names_to = "Dup_Type", values_to = "Proportion") %>% 
#       pivot_wider(names_from = minMax, values_from = Proportion)
#     
#     
#   } else if (i == "Low"){
#     tmp <- subset(sum_stats, CNV == i & is_nlr == "NLR")
#     low_ci <-
#       sapply(unique(tmp$Dup_Type),
#              function(k){
#                confint(lm(Count ~ 1, subset(tmp, Dup_Type == k)))
#              }, simplify = TRUE) %>% 
#       data.frame() %>%
#       setNames(nm = unique(tmp$Dup_Type)) %>% 
#       add_column(minMax = c("ymin", "ymax")) %>%
#       mutate(CNV = rep("Low", nrow(.))) %>%
#       pivot_longer(!c(minMax, CNV), names_to = "Dup_Type", values_to = "Proportion") %>% 
#       pivot_wider(names_from = minMax, values_from = Proportion)
#     
#   }
# }

############
### PLOT ###
############
# PLOT: Number of NLRs of each dup type for each genotype
nlr_sums <- 
  nlr %>%
  dplyr::select(Gene, Genotype, Dup_Type, Chr) %>% 
  group_by(Genotype, Dup_Type) %>% 
  tally()

pdf("figures/NLR_Duplication_highLow.pdf", height = 7, width = 7)
total_nlr %>%
  add_column(Dup_Type = rep("total", nrow(.)), .before = "Total") %>% 
  setNames(nm = c("Genotype", "Dup_Type", "n")) %>% 
  bind_rows(nlr_sums) %>% 
  distinct() %>% 
  add_column(CNV = ifelse(.$Genotype == "CCN-51" | .$Genotype == "ICS-1" |
                            .$Genotype == "GU-257E" | .$Genotype == "NA-246", "Low", "High")) %>%
  mutate(CNV = factor(CNV, 
                      levels = c("Low", "High"),
                      labels = c("Low", "High"))) %>%
  filter(Dup_Type != "singleton") %>%
  mutate(Dup_Type = factor(Dup_Type, 
                           levels = c("total", "dispersed", "proximal", "tandem", "segmental"),
                           labels = c("total", "dispersed", "proximal", "tandem", "segmental"))) %>%
  ggplot(aes(x = CNV, y = n, fill = CNV)) + 
  geom_jitter(width = 0.2, pch = 21, color = "black", size =4) + 
  facet_wrap(~Dup_Type, nrow = 3, ncol = 2) + 
  stat_summary(aes(color = CNV), fun = mean, shape = 5,
               position = position_dodge(width = 0.75),
               show.legend = FALSE, size = 1, color = "black") +
  stat_summary(aes(color = CNV),
               fun.data = mean_cl_normal,
               geom = "errorbar",
               width= 0.0,
               position = position_dodge(width = 0.75),
               show.legened = FALSE,
               color = "black") + 
  theme_minimal_grid(font_size = 15) + 
  scale_fill_manual(values = c("Low" = "#7876B1FF", "High" = "#E18727FF")) + 
  labs(x = NULL, y = "# NLR Genes", fill = NULL) + 
  theme(
    legend.position = c(0.85,0.15),
    legend.justification = "right",
    legend.direction = "vertical",
    strip.background = element_rect(color = "lightgrey")
  ) + 
  guides(fill=guide_legend(override.aes = list(size=8, shape = 21)))
dev.off()

# Model
mod <-
  nlr_sums %>%
  add_column(CNV = ifelse(.$Genotype == "CCN-51" | .$Genotype == "ICS-1" |
                            .$Genotype == "GU-257E" | .$Genotype == "NA-246", "Low", "High")) %>%
  filter(Dup_Type != "singleton")

# Fit a GLM with a negative binomial distribution--poisson didn't seem to fit
glm_mod <- glm.nb(n ~ Dup_Type + CNV + Dup_Type:CNV,  data = mod)

# Does a negative binomial dist fit the data? If p > 0.05 a negative binomial dist is appropriate
1 - pchisq(summary(glm_mod)$deviance, summary(glm_mod)$df.residual)

# Extrract pairwise comprisons using emmeans
mod <- summary(pairs(emmeans(glm_mod,  ~ Dup_Type + CNV + Dup_Type:CNV)))

mod <- 
  mod %>% 
  .[grepl("High.* - Low", .$contrast), ]

############
### Plot ### 
############

pdf("figures/NLR_DuplicationTypes_wholeGenome.pdf", width = 4.5)
x <- bind_rows(nlr, non_nlr) %>%
  group_by(is_nlr, Genotype, Dup_Type) %>%
  tally(name = "Count") %>%
  group_by(is_nlr, Genotype) %>%
  mutate(Proportion = Count/sum(Count)) %>%
  mutate(Dup_Type = factor(Dup_Type, levels = c("singleton","dispersed", "proximal", "tandem", "segmental"),
                           labels = c("singleton","dispersed", "proximal", "tandem", "segmental")))
  ggplot(aes(x = is_nlr, y = Proportion, fill = is_nlr))+ 
  geom_jitter( color = "black", pch = 21, 
               width = 0.05, size = 3, alpha = 1) +
  stat_summary(fun = mean, shape = 5, size = 1, show.legend = FALSE) +
  stat_summary(fun.data = mean_cl_normal,
               geom = "errorbar",
               width= 0.0,
               position = position_dodge(width = 0.75),
               show.legened = FALSE, 
               color = "black") +
  theme_minimal_grid(font_size = 15) + 
  labs(x = NULL, y = "% Relative to Whole", fill = NULL)+ 
  #scale_color_manual(values = c("#20b12d", "#5b91ff")) + 
  scale_y_continuous(labels=scales::percent_format(accuracy = 5L)) +
  scale_fill_manual(values = c("NLR" = "black", "Non-NLR" = "darkgrey")) +
  theme(
    legend.position = c(0.6, 0.15),
    legend.direction = "vertical",
    legend.text = element_text(size = 12), 
    axis.text.x = element_text(size = 15, angle = 45,  hjust = 1),
    axis.text.y = element_text(size = 15), 
    axis.title.y = element_text(size = 15), 
    strip.background = element_rect(color = "lightgrey")
  ) + 
  facet_wrap(~Dup_Type, ncol = 2, nrow = 3) + 
  guides(color = guide_legend(override.aes = list(size = 10)))
dev.off()

# Are the differences between groups, for each duplicate type, significantly different from one another? 
dup_sums <- 
  bind_rows(nlr, non_nlr) %>%
  group_by(is_nlr, Genotype, Dup_Type) %>%
  tally(name = "Count") %>%
  group_by(is_nlr, Genotype) %>%
  mutate(Proportion = Count/sum(Count)) 

mod <- aov(Proportion ~ is_nlr + Dup_Type + is_nlr:Dup_Type, dup_sums)
check_model(mod)
summary(mod)
TukeyHSD(mod)


############
### Plot ### 
############ 
# Duplicate Types Per Genotype  
sum_stats <- 
  nlr %>%
  group_by(Genotype, Dup_Type, Class) %>%
  tally(name = "Count") %>%
  group_by(Genotype) %>%
  mutate(Proportion = Count/sum(Count)) %>%
  filter(Class != "LRR") %>%
  mutate(Class = factor(Class,
                        levels = c("NL", "CNL", "TNL", "RNL"),
                        labels = c("NL", "CNL", "TNL", "RNL"))) %>%
  mutate(Dup_Type = factor(Dup_Type, 
                           levels = c("singleton", "dispersed", "proximal", "tandem", "segmental"),
                           labels = c("singleton", "dispersed", "proximal", "tandem", "segmental"))) %>%
  add_column(CNV = ifelse(.$Genotype == "CCN-51" | .$Genotype == "ICS-1" |
                            .$Genotype == "GU-257E" | .$Genotype == "NA-246", "Low", "High"))


for (i in unique(sum_stats$Class)){ 
  if (i == "NL"){
    tmp <- subset(sum_stats, Class == i)
    nl_ci <-
      sapply(unique(tmp$Dup_Type),
             function(k){
               confint(lm(Count ~ 1, subset(tmp, Dup_Type == k)))
             }, simplify = TRUE) %>% 
      data.frame() %>%
      setNames(nm = unique(tmp$Dup_Type)) %>% 
      add_column(minMax = c("ymin", "ymax")) %>%
      mutate(Class = rep(i, nrow(.))) %>% 
      pivot_longer(!c(minMax, Class), names_to = "Dup_Type", values_to = "Proportion") %>% 
      pivot_wider(names_from = minMax, values_from = Proportion)
    
    nl_ci[is.na(nl_ci)] <- 0
    }
  else if (i == "CNL"){
    tmp <- subset(sum_stats, Class == i)
    cnl_ci <-
      sapply(unique(tmp$Dup_Type),
             function(k){
               confint(lm(Count ~ 1, subset(tmp, Dup_Type == k)))
             }, simplify = TRUE) %>% 
      data.frame() %>%
      setNames(nm = unique(tmp$Dup_Type)) %>% 
      add_column(minMax = c("ymin", "ymax")) %>%
      mutate(Class = rep(i, nrow(.))) %>% 
      pivot_longer(!c(minMax, Class), names_to = "Dup_Type", values_to = "Proportion") %>% 
      pivot_wider(names_from = minMax, values_from = Proportion)
    
    cnl_ci[is.na(cnl_ci)] <- 0
    }
  else if (i == "TNL"){
    tmp <- subset(sum_stats, Class == i)
    tnl_ci <-
      sapply(unique(tmp$Dup_Type),
             function(k){
               confint(lm(Count ~ 1, subset(tmp, Dup_Type == k)))
             }, simplify = TRUE) %>% 
      data.frame() %>%
      setNames(nm = unique(tmp$Dup_Type)) %>% 
      add_column(minMax = c("ymin", "ymax")) %>%
      mutate(Class = rep(i, nrow(.))) %>% 
      pivot_longer(!c(minMax, Class), names_to = "Dup_Type", values_to = "Proportion") %>% 
      pivot_wider(names_from = minMax, values_from = Proportion)
    
    tnl_ci[is.na(tnl_ci)] <- 0
    }
  else if (i == "RNL"){
    tmp <- subset(sum_stats, Class == i)
    rnl_ci <-
      sapply(unique(tmp$Dup_Type),
             function(k){
               confint(lm(Count ~ 1, subset(tmp, Dup_Type == k)))
             }, simplify = TRUE) %>% 
      data.frame() %>%
      setNames(nm = unique(tmp$Dup_Type)) %>% 
      add_column(minMax = c("ymin", "ymax")) %>%
      mutate(Class = rep(i, nrow(.))) %>% 
      pivot_longer(!c(minMax, Class), names_to = "Dup_Type", values_to = "Proportion") %>% 
      pivot_wider(names_from = minMax, values_from = Proportion)
    
    rnl_ci[is.na(rnl_ci)] <- 0
  }
}

ci <- rbind(nl_ci, cnl_ci, tnl_ci, rnl_ci) %>%
  filter(Dup_Type != "singleton")

# Duplicate Types Per Genotype 
pdf("figures/NLR_DuplicationTypes_byClass.pdf", height = 4.5, width = 4.5)
nlr %>%
  group_by(Genotype, Dup_Type, Class) %>%
  tally(name = "Count") %>%
  group_by(Genotype) %>%
  mutate(Proportion = Count/sum(Count)) %>%
  filter(Class != "LRR") %>%
  filter(Dup_Type != "singleton") %>%
  mutate(Class = factor(Class,
                        levels = c("NL", "CNL", "TNL", "RNL"),
                        labels = c("NL", "CNL", "TNL", "RNL"))) %>%
  mutate(Dup_Type = factor(Dup_Type, 
                           levels = c("dispersed", "proximal", "tandem", "segmental"),
                           labels = c("dispersed", "proximal", "tandem", "segmental"))) %>%
  ggplot(aes(fill = Class, x = Class, y = Count)) +
  #geom_boxplot() +
  geom_jitter(width = 0.1, show.legend = TRUE, 
              alpha = 1, color = "black", pch = 21, 
              size = 4) +
  stat_summary(fun = mean, shape = 5, show.legend = FALSE, size = 1) +
  scale_fill_manual(values = c("#046C9A", "#D69C4E", "#ABDDDE", "#000000")) +
  theme_minimal_grid() + 
  labs(x = NULL, y = "# NLR Genes", fill = NULL) + 
  theme(
    legend.position = 'bottom',
    legend.justification = 'center',
    legend.direction = "horizontal",
    legend.text = element_text(size = 12), 
    axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 15), 
    axis.title.y = element_text(size = 15),
    strip.background = element_rect(color = "lightgrey")
  ) +
  stat_summary(fun.data = mean_cl_normal,
               geom = "errorbar",
               width= 0.0,
               position = position_dodge(width = 0.75),
               show.legened = FALSE, 
               color = "black") +
  facet_wrap(~Dup_Type, ncol = 2, nrow = 3) +
  guides(color=guide_legend(override.aes = list(size=4)))
dev.off()

# Model 
mod_df <-
  nlr %>%
  group_by(Genotype, Dup_Type, Class) %>%
  tally(name = "Count") %>%
  group_by(Genotype) %>%
  mutate(Proportion = Count/sum(Count)) %>%
  filter(Class != "LRR") %>%
  filter(Dup_Type != "singleton") %>%
  mutate(Class = factor(Class,
                        levels = c("NL", "CNL", "TNL", "RNL"),
                        labels = c("NL", "CNL", "TNL", "RNL"))) %>%
  mutate(Dup_Type = factor(Dup_Type, 
                           levels = c("dispersed", "proximal", "tandem", "segmental"),
                           labels = c("dispersed", "proximal", "tandem", "segmental"))) %>%
  add_column(CNV = ifelse(.$Genotype == "CCN-51" | .$Genotype == "ICS-1" |
                            .$Genotype == "GU-257E" | .$Genotype == "NA-246", "Low", "High"))

t.test(Count ~ CNV, mod_df)

############
### Plot ### 
############
ortho_sums <-
  ortho %>% 
  group_by(Orthogroup, Genotype) %>% 
  tally(name = "Count") %>% 
  add_column(CNV = ifelse(.$Genotype == "CCN-51" | .$Genotype == "ICS-1" |
                            .$Genotype == "GU-257E" | .$Genotype == "NA-246", "Low", "High")) %>%
  mutate(Orthogroup = as.character(Orthogroup)) %>%
  filter(Orthogroup == "4" | 
           Orthogroup == "16" | 
           Orthogroup == "6" | 
           Orthogroup == "567") %>%
  mutate(Orthogroup = paste0("OG", Orthogroup)) 

for (i in unique(ortho_sums$CNV)){ 
  if (i == "High"){
    tmp <- subset(ortho_sums, CNV == i)
    high_ci <-
      sapply(unique(tmp$Orthogroup),
             function(k){
               confint(lm(Count ~ 1, subset(tmp, Orthogroup == k)))
             }, simplify = TRUE) %>% 
      data.frame() %>% 
      add_column(minMax = c("ymin", "ymax")) %>%
      mutate(CNV = rep(i, nrow(.))) %>% 
      pivot_longer(!c(minMax, CNV), names_to = "Orthogroup", values_to = "Count") %>% 
      pivot_wider(names_from = minMax, values_from = Count)
  
  }
  else if (i == "Low"){
    tmp <- subset(ortho_sums, CNV == i)
    low_ci <-
      sapply(unique(tmp$Orthogroup),
             function(k){
               confint(lm(Count ~ 1, subset(tmp, Orthogroup == k)))
             }, simplify = TRUE) %>% 
      data.frame() %>% 
      add_column(minMax = c("ymin", "ymax")) %>%
      mutate(CNV = rep(i, nrow(.))) %>% 
      pivot_longer(!c(minMax, CNV), names_to = "Orthogroup", values_to = "Count") %>% 
      pivot_wider(names_from = minMax, values_from = Count)
    
  }
}
  
ci <- rbind(high_ci, low_ci)


# Plot the distribution of orthogroups
pdf("figures/NLR_Duplication_orthoCNV.pdf", width = 4.5, height = 4.5)
ortho %>% 
  group_by(Genotype, Orthogroup) %>% 
  tally(name = "Count") %>% 
  mutate(Orthogroup = as.character(Orthogroup)) %>%
  filter(Orthogroup == "4" | 
           Orthogroup == "16" | 
           Orthogroup == "6" | 
           Orthogroup == "567") %>%
  mutate(Orthogroup = paste0("OG", Orthogroup)) %>%
  add_column(CNV = ifelse(.$Genotype == "CCN-51" | .$Genotype == "ICS-1" |
                            .$Genotype == "GU-257E" | .$Genotype == "NA-246", "Low", "High")) %>%
  ggplot(aes(x = CNV, y = Count, color = CNV)) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2), 
             alpha = 0.6,
             size = 2) + 
  stat_summary(fun = mean, 
               shape = 5, 
               position = position_dodge(width = 0.75), 
               show.legend = FALSE, 
               size = 1) +
  scale_color_manual(values = c("black", "darkgrey")) + 
  labs(color = NULL, x = NULL, y = "# of NLR Genes") +
  theme_minimal_grid() +
  # geom_errorbar(data = ci, 
  #               aes(x = CNV, ymin = ymin, ymax = ymax, color = CNV),
  #               inherit.aes = FALSE, 
  #               position = position_dodge(width = 0.75), 
  #               width = 0.0,
  #               show.legend = FALSE) +
  stat_summary(aes(color = CNV), 
               fun.data = mean_cl_normal,  
               geom = "errorbar", 
               width= 0.0, 
               position = position_dodge(width = 0.75)) +
  theme(
    legend.position = 'bottom',
    legend.justification = 'center',
    legend.direction = "horizontal",
    legend.text = element_text(size = 12), 
    axis.text.x = element_text(size = 15, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 15, face = "bold"), 
    axis.title.y = element_text(size = 15, face = "bold"),
    strip.background = element_rect(color = "lightgrey")
  ) + 
  facet_wrap(~Orthogroup, scales = "free_y") + 
  guides(color=guide_legend(override.aes = list(size=4)))
dev.off()

# Model 
mod_df <- 
  ortho %>% 
  group_by(Genotype, Orthogroup) %>% 
  tally(name = "Count") %>% 
  mutate(Orthogroup = as.character(Orthogroup)) %>%
  filter(Orthogroup == "4" | 
           Orthogroup == "16" | 
           Orthogroup == "6" | 
           Orthogroup == "567") %>%
  mutate(Orthogroup = paste0("OG", Orthogroup)) %>%
  add_column(CNV = ifelse(.$Genotype == "CCN-51" | .$Genotype == "ICS-1" |
                            .$Genotype == "GU-257E" | .$Genotype == "NA-246", "Low", "High"))

mod <- aov(Count ~ Orthogroup + CNV + Orthogroup:CNV, mod_df)
#check_model(mod)
summary(mod)
mod <- TukeyHSD(mod)
mod <- data.frame(mod$`Orthogroup:CNV`)

############
### Plot ### 
############
# Plot average orthogroup size

pdf("figures/NLR_Duplication_orthoSize.pdf", width = 11)
ortho %>% 
  group_by(Orthogroup, Genotype) %>% 
  tally(name = "Count") %>%
  summarise(mCount = round(mean(Count), 1)) %>% 
  arrange(desc(mCount)) %>%
  mutate(Orthogroup = paste0("OG", Orthogroup)) %>% 
  mutate(Color = ifelse(Orthogroup == "OG4" | Orthogroup == "OG16" | 
                          Orthogroup == "OG6" | Orthogroup == "OG567", TRUE, FALSE)) %>%
  mutate(Lbl = ifelse(Color == TRUE, .$mCount, NA)) %>%
  #mutate(xAxis = ifelse(Color == TRUE, .$Orthogroup, "")) %>%
  ggplot(aes(x = reorder(Orthogroup, mCount), y = mCount, fill = Color)) + 
  geom_bar(stat = "identity", show.legend = FALSE) + 
  labs(x = NULL, y = "Average # of NLR Genes per Genotype") +
  theme_minimal_grid() +
  theme(
    axis.text.y = element_text(size = 20, face = "bold", hjust = 1,
                               color = c(rep("transparent", 23), rep("black", 4))),
    axis.text.x = element_text(size = 20, face = "bold"), 
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 20, face = "bold"),
    strip.background = element_rect(color = "lightgrey")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) + 
  scale_fill_manual(values = c("darkgrey", "#343c80")) + 
  geom_text(aes(label = Lbl), hjust = -0.05, size = 8) + 
  coord_flip()
dev.off()

############
### PLOT ###
############
# Models for orthogroup:chromosome 
gene_pos <- read.table(gene_pos_file, header = FALSE) %>%
  setNames(nm = c("Gene", "Chr", "Start", "Stop"))
ortho_pos <-
  ortho %>% 
  inner_join(gene_pos, by = "Gene") %>%
  group_by(Orthogroup, Genotype , Chr) %>% 
  tally(name = "Count") %>% 
  add_column(CNV = ifelse(.$Genotype == "CCN-51" | .$Genotype == "ICS-1" |
                            .$Genotype == "GU-257E" | .$Genotype == "NA-246", "Low", "High")) %>%
  mutate(Orthogroup = as.character(Orthogroup)) %>%
  filter(Orthogroup == "4" | 
           Orthogroup == "16" | 
           Orthogroup == "6" | 
           Orthogroup == "567") %>%
  mutate(Orthogroup = paste0("OG", Orthogroup)) %>% 
  mutate(Chr = gsub("^[A-z].*_Chr", "", Chr)) %>%
  mutate(Chr = gsub("chr", "", Chr)) %>% 
  filter(Chr != "Un_random") %>%
  .[-grep("Matina_scaffold.*", .$Chr),] %>%
  mutate(Chr = paste0("Chr", Chr)) %>%
  filter(Chr != "Chr0") %>%
  mutate(Chr = factor(Chr, 
                      levels = c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5",
                                      "Chr6", "Chr7", "Chr8", "Chr9", "Chr10"),
                      labels = c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5",
                                 "Chr6", "Chr7", "Chr8", "Chr9", "Chr10")))



mod <- lapply(split(ortho_pos, ortho_pos$Orthogroup), function(d) { aov(Count ~ Chr, data=d) })
lapply(mod, summary)
lapply(mod, function(k){
  data.frame(TukeyHSD(k)$Chr) %>% 
    dplyr::select(p.adj) %>% 
    rownames_to_column(var = "Comp") %>% 
    deframe() %>%
    multcompLetters(.)
})

############
### PLOT ###
############
pdf("figures/NLR_Duplication_orthoCount_pos.pdf")
ggplot(ortho_pos, aes(x = Chr, y = Count, color = CNV)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2), 
             alpha = 0.6, show.legend = FALSE) +
  stat_summary(aes(color = CNV), fun = mean, shape = 5,
               position = position_dodge(width = 0.75),
               show.legend = FALSE) +
  stat_summary(aes(color = CNV),
               fun.data = mean_cl_normal,
               geom = "errorbar",
               width= 0.0,
               position = position_dodge(width = 0.75),
               show.legened = FALSE) +
  scale_color_manual(values = c("black", "lightgrey")) +
  facet_wrap(~Orthogroup, scales = "free_x") + 
  theme_minimal_grid() +
  labs(x = NULL, y = "# NLR Genes", color = NULL) +
  theme(
    strip.background = element_rect(color = "lightgrey"), 
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom", 
    legend.direction = "horizontal", 
    legend.justification = "center"
  ) +
  guides(color=guide_legend(override.aes = list(size=4)))
dev.off()


############
### Plot ### 
############
pdf("~/Downloads/test.pdf")
inner_join(ortho, dup, by = "Gene") %>% 
  filter(Orthogroup == "4" | Orthogroup == "6" | 
           Orthogroup == "16" | Orthogroup == "567") %>%
  subset(Dup_Type != 'singleton') %>%
  group_by(Orthogroup, Genotype, Dup_Type) %>%
  tally(name = "Count") %>% 
  mutate(Orthogroup = paste0("OG", Orthogroup)) %>%
  add_column(CNV = ifelse(.$Genotype == "CCN-51" | .$Genotype == "ICS-1" |
                            .$Genotype == "GU-257E" | .$Genotype == "NA-246", "Low", "High")) %>%
  ggplot(aes(x = Dup_Type, y = Count, shape = CNV)) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2),
             show.legend = FALSE) + 
  stat_summary(aes(x = Dup_Type, y = Count), fun = mean,
               position = position_dodge(width = 0.75), 
               color = "red",
               show.legend = FALSE) +
  stat_summary(aes(x = Dup_Type, y = Count),
               fun.data = mean_cl_normal,
               geom = "errorbar",
               width= 0.0,
               position = position_dodge(width = 0.75),
               color = "red",
               show.legend = FALSE) +
  facet_wrap(~Orthogroup) + 
  theme_minimal_grid() +
  theme(
    axis.text.x = element_text(angle = 45, size = 15, face = "bold", hjust = 1),
    axis.text.y = element_text(size = 15, face = "bold", hjust = 1),
    axis.title.y = element_text(size = 15, face = "bold"),
    strip.background = element_rect(color = "lightgrey"), 
    legend.position = "bottom", 
    legend.direction = "horizontal", 
    legend.justification = "center"
  ) + 
  labs(x = NULL, y = "# of NLR Genes", color = NULL, alpha = NULL, shape = NULL) 
  #scale_color_manual(values = c("#046C9A", "#D69C4E", "#ABDDDE", "#000000")) + 
  #scale_alpha_manual(values = c(1, 0.4))
dev.off()

x <-
  inner_join(ortho, dup, by = "Gene") %>% 
  filter(Orthogroup == "4" | Orthogroup == "6" | 
           Orthogroup == "16" | Orthogroup == "567") %>%
  subset(Dup_Type != 'singleton') %>%
  group_by(Orthogroup, Genotype, Dup_Type) %>%
  tally(name = "Count") %>% 
  mutate(Orthogroup = paste0("OG", Orthogroup)) %>%
  add_column(CNV = ifelse(.$Genotype == "CCN-51" | .$Genotype == "ICS-1" |
                            .$Genotype == "GU-257E" | .$Genotype == "NA-246", "Low", "High"))

#mod <- lapply(split(x, x$Orthogroup), function(d) { aov(Count ~ CNV:Dup_Type, data=d) })
tmp <- subset(x, Orthogroup == "OG16")
mod <- aov(Count ~ CNV:Dup_Type, tmp)
x <- TukeyHSD(mod) %>%
  .[[1]] %>%
  data.frame(.) %>%
  dplyr::select(p.adj) %>%
  rownames_to_column(var = "Comp") %>%
  deframe()

multcompLetters(na.omit(x))

lapply(mod, summary)
lapply(mod, function(k){
  data.frame(TukeyHSD(k)$Chr) %>% 
    dplyr::select(p.adj) %>% 
    rownames_to_column(var = "Comp") %>% 
    deframe() %>%
    multcompLetters(.)
})

############
### Plot ### 
############
# What is the distribution of NLR duplication types across chromosomes? 
pdf("figures/NLR_DuplicationTypes_byChromosome.pdf", width = 12)
nlr %>%
  group_by(Genotype, Chr, Dup_Type) %>%
  tally(name = "Count") %>%
  mutate(Dup_Type = factor(Dup_Type, 
                           levels = c("singleton", "dispersed", "proximal", "tandem", "segmental"),
                           labels = c("singleton", "dispersed", "proximal", "tandem", "segmental"))) %>%
  mutate(Chr = factor(Chr, 
                          levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"),
                          labels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))) %>%
  ggplot(., aes(x = Dup_Type, y = Count, fill = Dup_Type)) + 
  geom_boxplot() +
  geom_point(position = position_dodge(width = 0.75)) +
  facet_wrap(~Chr) + 
  theme_minimal_grid() + 
  labs(x = NULL, y = "Number of NLRs", color = NULL, fill = "Duplication Type") + 
  theme(
    legend.position = "top",
    legend.text = element_text(size = 12), 
    axis.text.x = element_text(size = 15, face = "bold", angle = 90, vjust = 0.5, hjust = 0.75),
    axis.text.y = element_text(size = 20, face = "bold"), 
    axis.title.y = element_text(size = 20, face = "bold"),
    strip.background = element_rect(color = "lightgrey")
  )
dev.off()

############
### Plot ### 
############
# What is the distribution of NLR duplication types across genotype x chromosomes? 
pdf("figures/NLR_DuplicationTypes_byGenotypeChromosome.pdf", width = 50, height = 50)
nlr %>%
  group_by(Genotype, Chr, Dup_Type) %>%
  tally(name = "Count") %>%
  mutate(Dup_Type = factor(Dup_Type, 
                           levels = c("singleton", "dispersed", "proximal", "tandem", "segmental"),
                           labels = c("singleton", "dispersed", "proximal", "tandem", "segmental"))) %>%
  mutate(Chr = factor(Chr, 
                      levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"),
                      labels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))) %>%
  ggplot(., aes(x = Dup_Type, y = Count, fill = Dup_Type)) + 
  geom_boxplot() +
  geom_point(position = position_dodge(width = 0.75)) +
  facet_wrap(~Chr + Genotype) + 
  theme_minimal_grid() + 
  labs(x = NULL, y = "Number of NLRs", color = NULL, fill = "Duplication Type") + 
  theme(
    legend.position = "top",
    legend.text = element_text(size = 12), 
    axis.text.x = element_text(size = 15, face = "bold", angle = 90, vjust = 0.5, hjust = 0.75),
    axis.text.y = element_text(size = 20, face = "bold"), 
    axis.title.y = element_text(size = 20, face = "bold"),
    strip.background = element_rect(color = "lightgrey")
  )
dev.off()

############
### Plot ### 
############
# Plot NLR copy number variation vs duplication rates 
total_nlr <- 
  nlr %>%
  group_by(Genotype) %>% 
  tally(name = "Total")

p1 <- 
  nlr %>%
  group_by(Genotype, Dup_Type) %>% 
  tally(name = "Count") %>% 
  inner_join(total_nlr, by = "Genotype") %>% 
  subset(Dup_Type != "singleton") %>%
  mutate(Prop = Count / Total) %>% 
  ggplot(aes(x = Count, y = Total)) + 
  geom_point(color = "black") +
  theme_minimal_grid() +
  labs(x = "# of Duplicated Genes" , y = "# NLR Genes") + 
  theme(
    strip.background = element_rect(color = "lightgrey"),
    axis.text.x = element_text(angle = 45)
  ) + 
  geom_smooth(method = "lm", se = FALSE, linetype = "dotted") +
  #geom_label_repel(aes(x = Prop, y = Total, label = Genotype)) +
  facet_wrap(~Dup_Type, ncol = 4, nrow = 1, scales = "free")

p2 <- 
  nlr %>%
  group_by(Genotype, Dup_Type) %>% 
  tally(name = "Count") %>% 
  inner_join(total_nlr, by = "Genotype") %>% 
  subset(Dup_Type != "singleton") %>%
  mutate(Prop = Count / Total) %>% 
  ggplot(aes(x = Prop, y = Total)) + 
  geom_point(color = "black") +
  theme_minimal_grid() +
  labs(x = "Proportion of Duplicated Genes" , y = "# NLR Genes") + 
  theme(
    strip.background = element_rect(color = "lightgrey"),
    axis.text.x = element_text(angle = 45)
  ) + 
  geom_smooth(method = "lm", se = FALSE, linetype = "dotted") +
  #geom_label_repel(aes(x = Prop, y = Total, label = Genotype)) +
  facet_wrap(~Dup_Type, ncol = 4, nrow = 1, scales = "free")

pdf("figures/NLR_dupType_vs_CNV.pdf", height = 9)
plot_grid(p1,p2, ncol = 1, nrow = 2)
dev.off()

# Total gene content per genome 
gene_totals <- read.table(gene_content) %>%
  setNames(nm = c("Genotype", "Gene")) %>%
  inner_join(total_nlr, by = "Genotype") %>%
  setNames(nm = c("Genotype", "Total", "NLR")) %>%
  mutate(total_noNLR = Total - NLR) %>%
  dplyr::select(-Total) %>%
  pivot_longer(cols = c(NLR, total_noNLR), names_to = "Gene_Type", values_to = "Value") %>%
  mutate(Gene_Type = ifelse(.$Gene_Type == "NLR", "NLR", "Non-NLR"))

pdf("figures/total_gene_content_perGntype.pdf")
gene_totals %>% 
  mutate(Genotype = factor(Genotype, 
                           levels = c("CCN-51", "GU-257E", "ICS-1", "NA-246",
                                      "Matina", "Criollo", "IMC-105", "NA-807", "Pound-7",
                                      "SCA-6", "SPEC-54-1"),
                           
                           labels = c("CCN-51", "GU-257E", "ICS-1", "NA-246",
                                      "Matina", "Criollo", "IMC-105", "NA-807", "Pound-7",
                                      "SCA-6", "SPEC-54-1"))) %>% 
ggplot(aes(x = Genotype, y = Value, 
              fill = Gene_Type)) +
  geom_bar(stat = "identity", color = "black") + 
  theme_minimal_grid(font_size = 20) + 
  scale_fill_manual(values = c("#000000", "#FFFFFF")) + 
  scale_y_continuous(label = comma, limits=c(0, 30000), breaks = c(0, 10000, 20000, 30000)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    legend.justification = "center",
    legend.direction = "horizontal"
  ) + 
  labs(x = NULL, y = "Gene Number", fill = NULL) 
dev.off()

# Mod
gene_totals %>% 
  group_by(Genotype) %>% 
  summarise(num_genes = sum(Value)) %>%
  add_column(CNV = ifelse(.$Genotype == "CCN-51" | .$Genotype == "ICS-1" |
                            .$Genotype == "GU-257E" | .$Genotype == "NA-246", "Low", "High")) %>%
  wilcox.test(num_genes ~ CNV, .)

############
### PLOT ###
############
# AED scores for each genotype with highlighted NLRs
aed <- read.table(aed_scores, header = TRUE) %>%
  mutate(Genotype = gsub("_Chr[0-9].*", "", .$Gene)) %>%
  mutate(Gene = gsub("\\.1", "", Gene))

pdf("figures/aed_scores_perGntype.pdf")
ggplot(aed, aes(x = Genotype, y = AED)) + 
  geom_violin() +
  geom_boxplot(width = 0.1) + 
  theme_minimal_grid() +
  theme(
    axis.text.x = element_text(face = "bold", size = 15, angle = 45, hjust = 1),
    axis.text.y = element_text(face = "bold", size = 15),
    axis.title.y = element_text(face = "bold", size = 15),
    legend.position = c(0.3,0.98), 
    legend.direction = "horizontal"
  ) + 
  labs(x = NULL, y = "Gene Annotation Quality (AED)")
dev.off()


############
### PLOT ###
############
pdf("figures/aed_scores_nlrOnly_perGntype.pdf")
inner_join(aed, dplyr::select(nlr, -Genotype), by = "Gene") %>% 
  mutate(Genotype = factor(Genotype, 
                           levels = c("CCN-51", "GU-257E", "ICS-1", "NA-246",
                                      "Matina", "Criollo", "IMC-105", "NA-807", "Pound-7",
                                      "SCA-6", "SPEC-54-1"),
                           
                           labels = c("CCN-51", "GU-257E", "ICS-1", "NA-246",
                                      "Matina", "Criollo", "IMC-105", "NA-807", "Pound-7",
                                      "SCA-6", "SPEC-54-1"))) %>% 
  ggplot(aes(x = Genotype, y = AED)) + 
  geom_violin() +
  geom_boxplot(width = 0.1) + 
  theme_minimal_grid() +
  theme(
    axis.text.x = element_text(face = "bold", size = 15, angle = 45, hjust = 1),
    axis.text.y = element_text(face = "bold", size = 15),
    axis.title.y = element_text(face = "bold", size = 15),
    legend.position = c(0.3,0.98), 
    legend.direction = "horizontal"
  ) + 
  labs(x = NULL, y = "Gene Annotation Quality (AED)")
dev.off()

# Is _global_ AED dfferent between the high CNV and the low CNV varieties? 
inner_join(aed, dplyr::select(non_nlr, -Genotype), by = "Gene") %>%
  add_column(CNV = ifelse(.$Genotype == "CCN-51" | .$Genotype == "ICS-1" |
                            .$Genotype == "GU-257E" | .$Genotype == "NA-246", "Low", "High")) %>%
  t.test(AED ~ CNV, .)


inner_join(aed, dplyr::select(non_nlr, -Genotype), by = "Gene") %>%
  add_column(CNV = ifelse(.$Genotype == "CCN-51" | .$Genotype == "ICS-1" |
                            .$Genotype == "GU-257E" | .$Genotype == "NA-246", "Low", "High")) %>%
  group_by(CNV) %>% 
  summarise(mAED = mean(AED), 
            SD = sd(AED))

inner_join(aed, dplyr::select(non_nlr, -Genotype), by = "Gene") %>%
  add_column(CNV = ifelse(.$Genotype == "CCN-51" | .$Genotype == "ICS-1" |
                            .$Genotype == "GU-257E" | .$Genotype == "NA-246", "Low", "High")) %>% 
  ggplot(aes(x = AED, fill = CNV)) + 
  geom_density(alpha = 0.6) + 
  scale_fill_manual(values = c("steelblue", "firebrick")) + 
  theme_minimal_grid()

# Is AED different between the high CNV and low CNV varieties?
inner_join(aed, dplyr::select(nlr, -Genotype), by = "Gene") %>%
  add_column(CNV = ifelse(.$Genotype == "CCN-51" | .$Genotype == "ICS-1" |
                            .$Genotype == "GU-257E" | .$Genotype == "NA-246", "Low", "High")) %>%
  t.test(AED ~ CNV, .)

inner_join(aed, dplyr::select(nlr, -Genotype), by = "Gene") %>%
  add_column(CNV = ifelse(.$Genotype == "CCN-51" | .$Genotype == "ICS-1" |
                            .$Genotype == "GU-257E" | .$Genotype == "NA-246", "Low", "High")) %>%
  group_by(CNV) %>% 
  summarise(mAED = mean(AED),
            SD = sd(AED))

inner_join(aed, dplyr::select(nlr, -Genotype), by = "Gene") %>%
  add_column(CNV = ifelse(.$Genotype == "CCN-51" | .$Genotype == "ICS-1" |
                            .$Genotype == "GU-257E" | .$Genotype == "NA-246", "Low", "High")) %>% 
  ggplot(aes(x = AED, fill = CNV)) + 
  geom_density(alpha = 0.6) + 
  scale_fill_manual(values = c("steelblue", "firebrick")) + 
  theme_minimal_grid()
  

# Does total gene count differ between high CNV and low CNV varieties?
bind_rows(nlr, non_nlr) %>% 
  group_by(Genotype) %>%
  tally(name = "Count") %>% 
  add_column(CNV = ifelse(.$Genotype == "CCN-51" | .$Genotype == "ICS-1" |
                            .$Genotype == "GU-257E" | .$Genotype == "NA-246", "Low", "High")) %>% 
  aov(Count ~ CNV, .) %>% 
  summary()

mod <- aov(Count ~ CNV, x)
summary(mod)


###################
### PSEUDOGENES ###
###################
# Read in pseudogenes -- CCN-51, NA-246, and GU-257E all have 0 predicted pseudogenes
pseudos <- read.table("evol_analyses/pseudogene/allGntype_pseudogene_pos.txt",
                      header = TRUE) %>% 
  mutate(Genotype = gsub("_Chr.*", "", Chr)) %>% 
  mutate(Genotype = gsub("sCAffold_[0-9].*", "Matina", Genotype))

# Create data frames
pseudo_parents <- 
  pseudos %>% 
  separate(col = Gene, into = c("Pseudogene", "Gene"), sep = "\\_", extra = "merge")

dup_type <- read.table("evol_analyses/MCScanX/all_10x_pep_10aug2021_geneDupType.txt") %>% 
  setNames(nm = c("Gene", "dup_type")) %>% 
  mutate(dup_type = ifelse(dup_type == 0, "singleton",
                           ifelse(dup_type == 1, "dispersed",
                                  ifelse(dup_type == 2, "proximal",
                                         ifelse(dup_type == 3, "tandem", "segmental")))))

pseudo_parents <- inner_join(pseudo_parents, dup_type, by = "Gene") %>% 
  separate(col = Chr, into = c("Genotype", "Chr"), sep = "\\_", extra = "merge") %>% 
  mutate(Genotype = gsub("sCAffold", "Matina", Genotype))


############
### PLOT ###
############
pdf("figures/numPseudo_byCNV.pdf", width = 4, height = 4)
pseudo_parents %>% 
  group_by(Genotype) %>% 
  tally() %>% 
  rbind(., data.frame(
    Genotype = c("CCN-51", "NA-246", "GU-257E"),
    n = c(0, 0, 0)
  )) %>% 
  add_column(CNV = ifelse(.$Genotype == "CCN-51" | .$Genotype == "ICS-1" |
                           .$Genotype == "GU-257E" | .$Genotype == "NA-246", "Low", "High")) %>%
  mutate(CNV = factor(CNV, levels = c("Low", "High"), labels = c("Low", "High"))) %>%
  ggplot(aes(x = CNV, y = n, fill = CNV)) +
  scale_y_continuous(limits = c(0,200)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2), 
             alpha = 1, show.legend = TRUE, size = 4,
             pch = 21, color = "black") +
  stat_summary(color = "black", fun = mean, shape = 5,
               position = position_dodge(width = 0.75),
               show.legend = FALSE, 
               size = 1) +
  stat_summary(color = "black",
               fun.data = mean_cl_normal,
               geom = "errorbar",
               width= 0.0,
               position = position_dodge(width = 0.75),
               show.legened = FALSE) +
  scale_fill_manual(values = c("Low" = "#7876B1FF", "High" = "#E18727FF")) +
  theme_minimal_grid(font_size = 15) +
  labs(x = NULL, y = "# NLR Pseudogenes", fill = NULL) +
  theme(
    strip.background = element_rect(color = "lightgrey"), 
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom", 
    legend.direction = "horizontal", 
    legend.justification = "center"
  ) +
  guides(color=guide_legend(override.aes = list(size=10)))
dev.off()

# Mod 
pseudo_parents %>% 
  group_by(Genotype) %>% 
  tally() %>% 
  rbind(., data.frame(
    Genotype = c("CCN-51", "NA-246", "GU-257E"),
    n = c(0, 0, 0)
  )) %>%
  add_column(CNV = ifelse(.$Genotype == "CCN-51" | .$Genotype == "ICS-1" |
                                        .$Genotype == "GU-257E" | .$Genotype == "NA-246", "Low", "High")) %>% 
  wilcox.test(n ~ CNV, .)

############
### PLOT ###
############
pdf("figures/Pseudogene_dupType_prop.pdf", height = 4.5, width = 4.5)
pseudo_parents %>%
  group_by(Genotype) %>% 
  add_tally(name = "total") %>% 
  group_by(Genotype, dup_type) %>% 
  add_tally() %>% 
  mutate(proportion = n / total) %>% 
  dplyr::select(Genotype, dup_type, proportion) %>% 
  distinct() %>% 
  add_column(CNV = ifelse(.$Genotype == "CCN-51" | .$Genotype == "ICS-1" |
                            .$Genotype == "GU-257E" | .$Genotype == "NA-246", "Low", "High")) %>% 
  ggplot(., aes(x = dup_type, y = proportion, color =CNV)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2), 
             alpha = 0.6, show.legend = FALSE, size = 3) +
  stat_summary(aes(color = CNV), fun = mean, shape = 5,
               position = position_dodge(width = 0.75),
               show.legend = FALSE, 
               size = 1) +
  stat_summary(aes(color = CNV),
               fun.data = mean_cl_normal,
               geom = "errorbar",
               width= 0.0,
               position = position_dodge(width = 0.75),
               show.legened = FALSE) +
  scale_color_manual(values = c("black", "lightgrey")) +
  theme_minimal_grid() +
  labs(x = NULL, y = "Proportion of Total Pseudogenes", color = NULL) +
  theme(
    strip.background = element_rect(color = "lightgrey"), 
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom", 
    legend.direction = "horizontal", 
    legend.justification = "center"
  ) +
  guides(color=guide_legend(override.aes = list(size=4)))
dev.off()


############
### PLOT ###
############
## PLOT: NLR Density
pdf("figures/NLR_Duplication_pos.pdf", height = 4.5, width = 4.5)
nlr %>%
  group_by(Chr, Genotype) %>% 
  tally() %>%
  mutate(Chr = gsub("^", "Chr", Chr)) %>%
  mutate(Chr = factor(Chr,
                      labels = c("Chr0", "Chr1", "Chr2", "Chr3", "Chr4",
                                 "Chr5", "Chr6", "Chr7", "Chr8", "Chr9", "Chr10"),
                      levels = c("Chr0", "Chr1", "Chr2", "Chr3", "Chr4",
                                 "Chr5", "Chr6", "Chr7", "Chr8", "Chr9", "Chr10"))) %>%
  # mutate(Chr = factor(Chr,
  #                     labels = c("0", "1", "2", "3", "4",
  #                                "5", "6", "7", "8", "9", "10"),
  #                     levels = c("0", "1", "2", "3", "4",
  #                                "5", "6", "7", "8", "9", "10"))) %>%
  add_column(CNV = ifelse(.$Genotype == "CCN-51" | .$Genotype == "ICS-1" |
                            .$Genotype == "GU-257E" | .$Genotype == "NA-246", "Low", "High")) %>% 
  mutate(CNV = factor(CNV, levels = c("Low", "High"), labels = c("Low", "High"))) %>%
  ggplot(., aes(x = Chr, y = n, fill = CNV)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2), 
             alpha = 1, show.legend = FALSE, pch = 21,
             color = "black", size = 4) +
  stat_summary(color = "black", fun = mean, shape = 5,
               position = position_dodge(width = 0.75),
               show.legend = FALSE, size = 1) +
  stat_summary(color = "black",
               fun.data = mean_cl_normal,
               geom = "errorbar",
               width= 0.0,
               position = position_dodge(width = 0.75),
               show.legened = FALSE) +
  scale_fill_manual(values = c("Low" = "#7876B1FF", "High" = "#E18727FF")) +  
  theme_minimal_grid(font_size = 15) +
  labs(x = NULL, y = "# NLR Genes", color = NULL) +
  theme(
    strip.background = element_rect(color = "lightgrey"), 
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom", 
    legend.direction = "horizontal", 
    legend.justification = "center"
  ) +
  guides(color=guide_legend(override.aes = list(size=4)))
dev.off()


############
### PLOT ###
############
## PLOT: NLR Density x duplication type
pdf("figures/NLR_pos_dupType.pdf", height = 9, width = 9)
nlr %>%
  group_by(Chr, Genotype, Dup_Type) %>% 
  tally() %>%
  mutate(Chr = gsub("^", "Chr", Chr)) %>%
  mutate(Chr = factor(Chr,
                      labels = c("Chr0", "Chr1", "Chr2", "Chr3", "Chr4",
                                 "Chr5", "Chr6", "Chr7", "Chr8", "Chr9", "Chr10"),
                      levels = c("Chr0", "Chr1", "Chr2", "Chr3", "Chr4",
                                 "Chr5", "Chr6", "Chr7", "Chr8", "Chr9", "Chr10"))) %>%
  subset(Dup_Type != "singleton") %>%
  # mutate(Chr = factor(Chr,
  #                     labels = c("0", "1", "2", "3", "4",
  #                                "5", "6", "7", "8", "9", "10"),
  #                     levels = c("0", "1", "2", "3", "4",
  #                                "5", "6", "7", "8", "9", "10"))) %>%
  add_column(CNV = ifelse(.$Genotype == "CCN-51" | .$Genotype == "ICS-1" |
                            .$Genotype == "GU-257E" | .$Genotype == "NA-246", "Low", "High")) %>% 
  mutate(CNV = factor(CNV, levels = c("Low", "High"), labels = c("Low", "High"))) %>%
  ggplot(., aes(x = Chr, y = n, fill = CNV)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2), 
             alpha = 1, show.legend = FALSE, size = 3,
             pch = 21, color = "black") +
  stat_summary(aes(color = CNV), fun = mean, shape = 5,
               position = position_dodge(width = 0.75),
               show.legend = FALSE, size = 1, color = "black") +
  stat_summary(aes(color = CNV),
               fun.data = mean_cl_normal,
               geom = "errorbar",
               width= 0.0,
               position = position_dodge(width = 0.75),
               show.legened = FALSE,
               color = "black") +
  scale_fill_manual(values = c("Low" = "#7876B1FF", "High" = "#E18727FF")) +
  theme_minimal_grid(font_size = 20) +
  labs(x = NULL, y = "# NLR Genes", color = NULL) +
  theme(
    strip.background = element_rect(color = "lightgrey"), 
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom", 
    legend.direction = "horizontal", 
    legend.justification = "center"
  ) +
  guides(color=guide_legend(override.aes = list(size=4))) +
  facet_wrap(~Dup_Type)
dev.off()

# Where do we see the largest differences between low vs high, and what types of duplication explain this?
dupType_pos <- 
  nlr %>%
  group_by(Chr, Genotype, Dup_Type) %>% 
  tally() %>%
  mutate(Chr = gsub("^", "Chr", Chr)) %>%
  mutate(Chr = factor(Chr,
                      labels = c("Chr0", "Chr1", "Chr2", "Chr3", "Chr4",
                                 "Chr5", "Chr6", "Chr7", "Chr8", "Chr9", "Chr10"),
                      levels = c("Chr0", "Chr1", "Chr2", "Chr3", "Chr4",
                                 "Chr5", "Chr6", "Chr7", "Chr8", "Chr9", "Chr10"))) %>%
  subset(Dup_Type != "singleton") %>%
  # mutate(Chr = factor(Chr,
  #                     labels = c("0", "1", "2", "3", "4",
  #                                "5", "6", "7", "8", "9", "10"),
  #                     levels = c("0", "1", "2", "3", "4",
  #                                "5", "6", "7", "8", "9", "10"))) %>%
  add_column(CNV = ifelse(.$Genotype == "CCN-51" | .$Genotype == "ICS-1" |
                            .$Genotype == "GU-257E" | .$Genotype == "NA-246", "Low", "High")) %>% 
  mutate(CNV = factor(CNV, levels = c("Low", "High"), labels = c("Low", "High"))) %>%
  group_by(Chr, Dup_Type, CNV) %>%
  summarise(mGene = mean(n)) %>%
  pivot_wider(names_from = CNV, values_from = mGene) %>%
  mutate(High = ifelse(is.na(High), 0, High),
         Low = ifelse(is.na(Low), 0, Low)) %>% 
  mutate(effectSize = High - Low)


############
### PLOT ###
############
## PLOT: Pseudogene density
pdf("figures/NLRPseudogene_Duplication_pos.pdf", height = 4.5, width = 4.5)

tmp <- 
  pseudo_parents %>% 
  .[grep("Chr", .$Chr),] %>%
  group_by(Chr, Genotype) %>%
  tally() %>% 
  mutate(Chr = factor(Chr, 
                      levels = c("Chr0", "Chr1", "Chr2", "Chr3", "Chr4",
                                 "Chr5", "Chr6", "Chr7", "Chr8", "Chr9", "Chr10"),
                      labels = c("Chr0", "Chr1", "Chr2", "Chr3", "Chr4",
                                 "Chr5", "Chr6", "Chr7", "Chr8", "Chr9", "Chr10"))) %>%
  add_column(CNV = ifelse(.$Genotype == "CCN-51" | .$Genotype == "ICS-1" |
                            .$Genotype == "GU-257E" | .$Genotype == "NA-246", "Low", "High")) %>%
  mutate(CNV = factor(CNV, levels = c("Low", "High"), labels = c("Low", "High"))) %>% 
  bind_rows(., 
        data.frame(
          Chr = rep(paste0("Chr", 0:10), 3),
          Genotype = rep(c("CCN-51", "NA-246", "GU-257E"), 11),
          n = rep(c(0, 0, 0), 11),
          CNV = rep("Low", 33)
  ))

tmp %>% 
  ggplot(., aes(x = Chr, y = n, fill = CNV)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2), 
             alpha = 1, show.legend = TRUE, size = 4,
             pch = 21, color = "black") +
  stat_summary(color = "black", fun = mean, shape = 5,
               position = position_dodge(width = 0.75),
               show.legend = FALSE, size = 1) +
  stat_summary(color = "black",
               fun.data = mean_cl_normal,
               geom = "errorbar",
               width= 0.0,
               position = position_dodge(width = 0.75),
               show.legened = FALSE) +
  scale_fill_manual(values = c("Low" = "#7876B1FF", "High" = "#E18727FF")) +
  theme_minimal_grid() +
  labs(x = NULL, y = "# NLR Pseudogenes", fill = NULL) +
  theme(
    strip.background = element_rect(color = "lightgrey"), 
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom", 
    legend.direction = "horizontal", 
    legend.justification = "center"
  ) +
  scale_y_continuous(limits = c(0,50)) +
  guides(color=guide_legend(override.aes = list(size=10, shape = 21)))
dev.off()

###########
### TEs ###
###########

############
### PLOT ###
############
## PLOT: TE Density
te_cnts <- te %>% 
  separate(col = Chr, into = c("Genotype", "Chr"), sep = "\\_", extra = "merge") %>% 
  group_by(Chr, Genotype) %>% 
  tally()

pdf("figures/TE_density_pos.pdf", height = 4.5, width = 4.5)
te_cnts %>%
  mutate(Chr = factor(Chr, 
                      levels = c("Chr0", "Chr1", "Chr2", "Chr3", "Chr4",
                                 "Chr5", "Chr6", "Chr7", "Chr8", "Chr9", "Chr10"),
                      labels = c("Chr0", "Chr1", "Chr2", "Chr3", "Chr4",
                                 "Chr5", "Chr6", "Chr7", "Chr8", "Chr9", "Chr10"))) %>%
  ggplot(., aes(x = Chr, y = n)) +
  geom_jitter(width = 0.2, 
              alpha = 0.6, show.legend = FALSE) +
  stat_summary(color = "black", fun = mean, shape = 5,
               position = position_dodge(width = 0.75),
               show.legend = FALSE) +
  stat_summary(color = "black",
               fun.data = mean_cl_normal,
               geom = "errorbar",
               width= 0.0,
               position = position_dodge(width = 0.75),
               show.legened = FALSE) +
  #scale_color_manual(values = c("black", "lightgrey")) +
  theme_minimal_grid() +
  labs(x = NULL, y = "# LTRs", color = NULL) +
  theme(
    strip.background = element_rect(color = "lightgrey"), 
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom", 
    legend.direction = "horizontal", 
    legend.justification = "center"
  ) +
  guides(color=guide_legend(override.aes = list(size=4)))
dev.off()





