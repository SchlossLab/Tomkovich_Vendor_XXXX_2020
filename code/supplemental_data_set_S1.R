source("code/functions.R")

#Generate single excel workbook containing sheets for Data Set S1

#Sheets 1 and 2: Alpha diversity Kruskal-Wallis and pairwise wilcoxon results----
#Combine Kruskal-Wallis test results for shannon and richness diversity metrics
shannon_results <- read_tsv("data/process/shannon_stats_all_days.tsv") %>% 
  mutate(diversity_metric = "Shannon Diversity")
richness_results <- read_tsv("data/process/richness_stats_all_days.tsv") %>% 
  mutate(diversity_metric = "Richness")

#Sheet 1: Diversity kruskal-Wallis test results with adjusted P values
diversity_kruskal_stats_adjust <- rbind(shannon_results, richness_results)

#Combine pairwise wilcoxon test results for shannon and richness diversity metrics
shannon_pairwise_results <- read_tsv("data/process/shannon_stats_sig_days.tsv") %>% 
  mutate(diversity_metric = "Shannon Diversity")
richness_pairwise_results <- read_tsv("data/process/richness_stats_sig_days.tsv") %>% 
  mutate(diversity_metric = "Richness")
#Sheet 2: Diversity pairwise wilcoxon results
diversity_pairwise_results <- rbind(shannon_pairwise_results, richness_pairwise_results)

#Write PERMANOVA results to data set excel file.
#Sheet 3: PERMANOVA results for day -1, day 0, and day 1----
adonis_dn1_to_0 <- read_tsv("data/process/adonis_dn1-1.tsv")
#Sheet 4: PERMANOVA results within sources at day -1----
dn1_source <- read_tsv("data/process/adonis_dn1_source.tsv")

#OTU statistics for differences across sources of mice on day -1, 0 and 1----
otu_stats_dn1to1_combined <- read_tsv("data/process/otu_stats_dn1to1_combined.tsv")
#Sheet 5: OTUs that vary between sources on day -1----
kw_w_sig_dn1 <- otu_stats_dn1to1_combined %>% filter(day == -1)

#Sheet 6: C. difficile CFU statistical results----
cfu_kruskal_stats_adjust <- read_tsv("data/process/cfu_stats_all_days.tsv")

#Sheet 7: Mouse weight change statistical results----
weight_kruskal_stats_adjust <- read_tsv("data/process/weight_stats_all_days.tsv")

#Sheet 8: OTUs that vary between sources on day 0----
kw_w_sig_d0 <- otu_stats_dn1to1_combined %>% filter(day == 0)

#Sheet 9: OTUs that changed after clindamycin treatment----
o_dn1to0_pairs_stats_adjust <- read_tsv("data/process/otu_stats_dn1to0.tsv") 

#Sheet 10: OTUs that vary between sources on day 1----
kw_w_sig_d1 <- otu_stats_dn1to1_combined %>% filter(day == 1)

#Sheet 11: PERMANOVA results all sources over time
adonis_all_samples <- read_tsv("data/process/adonis_all.tsv")

#Sheets 12-13: Statistical results of L2-regularized logistic regression model performances----
#Sheet 12: Compare test AUCs for each model to 0.5. Wilcoxan signed rank test. Null hypothesis is the distribution of test AUCs for each model is symmetric about mu (0.5)
all_vs_random <-  read_tsv("data/process/classification_to_random.tsv") 
#Sheet 13: Pairwise comparisons of test AUCs from the 3 models
model_stats_pairwise <- read_tsv("data/process/classification_model_pairwise_stats.tsv") 

#Sheet 14: Top 20 most important OTUs for each of the 3 models----
combined_top20_all_models <- read_tsv("data/process/combined_top20_otus_all_models.tsv") 

#Sheet 15: OTUs that varied across sources over at least 1 timepoint----
otu_stats_dn1to9_combined <- read_tsv("data/process/otu_kw_stats_dn1to9.tsv") 

#Combine all the above tables as sheets---
sheets_1_15 <- list("Sheet 1" = diversity_kruskal_stats_adjust, 
                            "Sheet 2" = diversity_pairwise_results,
                            "Sheet 3" = adonis_dn1_to_0,
                            "Sheet 4"  = dn1_source, 
                            "Sheet 5" = kw_w_sig_dn1, 
                            "Sheet 6" = cfu_kruskal_stats_adjust,
                            "Sheet 7" = weight_kruskal_stats_adjust,
                            "Sheet 8" = kw_w_sig_d0, 
                            "Sheet 9" = o_dn1to0_pairs_stats_adjust,
                            "Sheet 10" = kw_w_sig_d1,
                            "Sheet 11" = adonis_all_samples,
                            "Sheet 12" = all_vs_random, 
                            "Sheet 13" = model_stats_pairwise,
                            "Sheet 14" = combined_top20_all_models,
                            "Sheet 15" = otu_stats_dn1to9_combined)

#Write sheets to  Data Set S1 excel file----
write_xlsx(sheets_1_15, "submission/Data Set S1.xlsx", format_headers = TRUE)
