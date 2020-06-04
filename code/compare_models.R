source("code/functions.R")

#Read in OTU level models
dayn1_file <- list.files(path= 'data/process/classification', pattern='combined_best_hp_results_dayn1_60', full.names = TRUE) 
otu_dayn1 <- map_df(dayn1_file, read_files) %>% 
  mutate(model_name = "otu_dayn1") 
day0_file <- list.files(path= 'data/process/classification', pattern='combined_best_hp_results_day0_60', full.names = TRUE) 
otu_day0 <- map_df(day0_file, read_files) %>% 
  mutate(model_name = "otu_day0") 
day1_file <- list.files(path= 'data/process/classification', pattern='combined_best_hp_results_day1_60', full.names = TRUE) 
otu_day1 <- map_df(day1_file, read_files) %>% 
  mutate(model_name = "otu_day1") 

#Read in family level models
dayn1_file <- list.files(path= 'data/process/classification', pattern='combined_best_hp_results_dayn1_60_family', full.names = TRUE) 
family_dayn1 <- map_df(dayn1_file, read_files) %>% 
  mutate(model_name = "family_dayn1") 
day0_file <- list.files(path= 'data/process/classification', pattern='combined_best_hp_results_day0_60_family', full.names = TRUE) 
family_day0 <- map_df(day0_file, read_files) %>% 
  mutate(model_name = "family_day0") 
day1_file <- list.files(path= 'data/process/classification', pattern='combined_best_hp_results_day1_60_family', full.names = TRUE) 
family_day1 <- map_df(day1_file, read_files) %>% 
  mutate(model_name = "family_day1")

#Combine the 6 dataframes for the 3 OTU level and 3 family level models 
all <- rbind(otu_dayn1, otu_day0, otu_day1, family_dayn1, family_day0, family_day1)

#Compare test AUCs for each model to 0.5. Wilcoxan signed rank test. Null hypothesis is the distribution of test AUCs for each model is symmetric about mu (0.5)
all %>% 
  group_by(model_name) %>% 
  summarize(test=list(tidy(wilcox.test(test_aucs, mu=0.5, paired = FALSE, alternative="greater"))), 
            summary=list(tidy(summary(test_aucs)))) %>% 
  unnest() %>% 
  write_tsv("data/process/classification_to_random.tsv") %>% 
  #Also write results to supplemental table excel file
  write_xlsx("submission/table_S12_classification_to_random.xlsx", format_headers = FALSE)
#All models perform significantly better than random (AUC of 0.5)

#Kruskal-wallis test comparing test AUCs from all 6 models: 
set.seed(19881117) #Same seed used for mothur analysis
kruskal.test(test_aucs~model_name, data = all) 
#Kruskal-Wallis rank sum test. p-value < 2.2e-16

#Pairwise comparisons of test AUCs from the 6 models
model_stats_pairwise <- 
  pairwise.wilcox.test(g = as.factor(all$model_name), x=all$test_aucs, p.adjust.method = "BH") %>% 
  tidy() %>% 
  mutate(compare=paste(group1, group2, sep="-")) %>%
  rename(p.value.adj = p.value) %>% #Specify that p.value is actually Benjamini-Hochberg adjusted P value
  select(compare, p.value.adj) %>% 
  rename(comparison = compare) %>% 
  arrange(p.value.adj) %>% 
  write_tsv("data/process/classification_model_pairwise_stats.tsv") %>% 
  #Also write results to supplemental table excel file
  write_xlsx("submission/table_S13_classification_model_pairwise_stats.xlsx", format_headers = FALSE)

#Compare each model's cv and test AUCs  
all %>% 
  group_by(model_name) %>% 
  nest() %>% 
  mutate(test = map(data, ~tidy(wilcox.test(x=.x$cv_aucs, y=.x$test_aucs, exact=FALSE))),
         cv_lci = map(data, ~quantile(.x$cv_aucs, 0.025)),
         cv_q1 = map(data, ~quantile(.x$cv_aucs, 0.25)),
         cv_median = map(data, ~median(.x$cv_aucs)),
         cv_q3 = map(data, ~quantile(.x$cv_aucs, 0.75)),
         cv_uci = map(data, ~quantile(.x$cv_aucs, 0.975)),
         test_lci = map(data, ~quantile(.x$test_aucs, 0.025)),
         test_q1 = map(data, ~quantile(.x$test_aucs, 0.25)),
         test_median = map(data, ~median(.x$test_aucs)),
         test_q3 = map(data, ~quantile(.x$test_aucs, 0.75)),
         test_uci = map(data, ~quantile(.x$test_aucs, 0.975))
  ) %>% 
  select(model_name, test,
         cv_lci, cv_q1, cv_median, cv_q3, cv_uci,
         test_lci, test_q1, test_median, test_q3, test_uci
  ) %>%
  unnest() %>%
  write_tsv("data/process/classification_cv_test_compare.tsv")
