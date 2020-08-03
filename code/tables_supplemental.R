source("code/functions.R")

#Generate single excel workbook containing sheets of all suppplemental tables

#Table S1 and S2: Alpha diversity Kruskal-Wallis and pairwise wilcoxon results----
#Combine Kruskal-Wallis test results for shannon and richness diversity metrics
shannon_results <- read_tsv("data/process/shannon_stats_all_days.tsv") %>% 
  mutate(diversity_metric = "Shannon Diversity")
richness_results <- read_tsv("data/process/richness_stats_all_days.tsv") %>% 
  mutate(diversity_metric = "Richness")

#Table S1: Diversity kruskal-Wallis test results with adjusted P values
diversity_kruskal_stats_adjust <- rbind(shannon_results, richness_results)

#Combine pairwise wilcoxon test results for shannon and richness diversity metrics
shannon_pairwise_results <- read_tsv("data/process/shannon_stats_sig_days.tsv") %>% 
  mutate(diversity_metric = "Shannon Diversity")
richness_pairwise_results <- read_tsv("data/process/richness_stats_sig_days.tsv") %>% 
  mutate(diversity_metric = "Richness")
#Table S2: Diversity pairwise wilcoxon results
diversity_pairwise_results <- rbind(shannon_pairwise_results, richness_pairwise_results)

#Write PERMANOVA results to supplemental table excel file.
#Table S3: PERMANOVA results for day -1, day 0, and day 1----
adonis_dn1_to_0 <- read_tsv("data/process/adonis_dn1-1.tsv")
#Table S4: PERMANOVA results within sources at day -1----
dn1_source <- read_tsv("data/process/adonis_dn1_source.tsv")

#OTU statistics for differences across sources of mice on day -1, 0 and 1----
otu_stats_dn1to1_combined <- read_tsv("data/process/otu_stats_dn1to1_combined.tsv")
#Table S5: OTUs that vary between sources on day -1----
kw_w_sig_dn1 <- otu_stats_dn1to1_combined %>% filter(day == -1)

#Table S6: C. difficile CFU statistical results----
cfu_kruskal_stats_adjust <- read_tsv("data/process/cfu_stats_all_days.tsv")

#Table S7: Mouse weight change statistical results----
weight_kruskal_stats_adjust <- read_tsv("data/process/weight_stats_all_days.tsv")

#Table S8: OTUs that vary between sources on day 0----
kw_w_sig_d0 <- otu_stats_dn1to1_combined %>% filter(day == 0)

#Table S9: OTUs that changed after clindamycin treatment----
o_dn1to0_pairs_stats_adjust <- read_tsv(path = "data/process/otu_stats_dn1to0.tsv") 

#Table S10: OTUs that vary between sources on day 1----
kw_w_sig_d1 <- otu_stats_dn1to1_combined %>% filter(day == 1)

#Table S11: PERMANOVA results all sources over time
adonis_all_samples <- read_tsv("data/process/adonis_all.tsv")

#Table S12-S13: Statistical results of L2-regularized logistic regression model performances----
#Table S12: Compare test AUCs for each model to 0.5. Wilcoxan signed rank test. Null hypothesis is the distribution of test AUCs for each model is symmetric about mu (0.5)
all_vs_random <-  read_tsv("data/process/classification_to_random.tsv") 
#Table S13: Pairwise comparisons of test AUCs from the 3 models
model_stats_pairwise <- read_tsv("data/process/classification_model_pairwise_stats.tsv") 

#Table S14: Top 20 most important OTUs for each of the 3 models----
combined_top20_all_models <- read_tsv("data/process/combined_top20_otus_all_models.tsv") 

#Table S15: OTUs that varied across sources over at least 1 timepoint----
otu_stats_dn1to9_combined <- read_tsv("data/process/otu_kw_stats_dn1to9.tsv") 

#Combine all the above tables as sheets---
table_S1_S15_sheets <- list("S1_Diversity_Kruskal_Wallis" = diversity_kruskal_stats_adjust, "S2_diversity_pairwise" = diversity_pairwise_results,
                            "S3_PERMANOVAdn1_d0_d1" = adonis_dn1_to_0, "S4_PERMANOVA_source_dn1"  = dn1_source, 
                            "S5_source_OTUs_day_-1" = kw_w_sig_dn1, 
                            "S6_CFU_stats" = cfu_kruskal_stats_adjust, "S7_weight_stats" = weight_kruskal_stats_adjust,
                            "S8_source_OTUs_day_0" = kw_w_sig_d0, "S9_clindamycin_OTUs" = o_dn1to0_pairs_stats_adjust,
                            "S10_source_OTUs_day_1" = kw_w_sig_d1,
                            "S11_PERMANOVA_all_timepoints" = adonis_all_samples,
                            "S12_models_vs_random" = all_vs_random, "S13_models_pairwise" = model_stats_pairwise,
                            "S14_top_OTUs_per_model" = combined_top20_all_models,
                            "S15_source_OTUs_over_time" = otu_stats_dn1to9_combined)

#Write supplemental table sheets to excel file----
write_xlsx(table_S1_S15_sheets, "submission/tables_S1-S15.xlsx", format_headers = TRUE)
