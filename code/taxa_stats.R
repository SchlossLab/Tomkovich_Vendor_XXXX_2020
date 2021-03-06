source("code/functions.R")

# Import taxonomy into data frame and clean up taxonomy names
taxonomy <- read_tsv(file="data/process/vendors.taxonomy") %>%
  rename_all(tolower) %>% #remove uppercase from column names
  # Split taxonomic information into separate columns for each taxonomic level  
  mutate(taxonomy=str_replace_all(taxonomy, c("\\(\\d*\\)" = "", #drop digits with parentheses around them
                                              ';$' = "", #removes semi-colon at end of line
                                              'Bacteria_unclassified' = 'Unclassified',
                                              "Clostridium_" = "Clostridium ", #Remove underscores after Clostridium
                                              "_" = " ", #Removes all other underscores
                                              "unclassified" = "Unclassified"))) %>% 
  # Separate taxonomic levels into separate columns according to semi-colon.
  separate(taxonomy, into=c("kingdom", "phylum", "class", "order", "family", "genus"), sep=';')

# Import otu_data for samples
otu_data <- read_tsv("data/process/vendors.subsample.shared", col_types=cols(Group=col_character())) %>% 
  select(-label, -numOtus) %>% 
  rename(id=Group) %>% 
  gather(-id, key="otu", value="count") %>% 
  filter(!otu== "Otu0020") %>% #remove C. difficile (Otu0020) from the input data since C. difficile colonization status at day 7 is more accurately quantified via plating 
  mutate(rel_abund=count/5437) #Use 5437, because this is the subsampling parameter chosen.)

#Merge otu_data to taxonomy data frame
agg_taxa_data <- inner_join(otu_data, taxonomy)

# Function to summarize relative abundance level for a given taxonomic level (ex. genus, family, phlyum, etc.)
agg_taxonomic_data <- function(taxonomic_level) {
  agg_taxa_data %>% 
    group_by(id, {{ taxonomic_level }}) %>% #Embracing treats the taxonomic_level argument as a column name
    summarize(agg_rel_abund=sum(rel_abund)) %>% 
    # Merge relative abundance data to specifci taxonomic_level data
    inner_join(., metadata, by = "id") %>% 
    ungroup() 
}

# Relative abundance data at the otu level:
agg_otu_data <- agg_taxonomic_data(otu)

# Number of detected OTUs in the dataset:
detected_otus <- agg_otu_data %>% 
  distinct(otu) %>% count #Selects the number of unique OTUs detected in dataset

#Rename otus to match naming convention used in class_interpretation.R (interpretation of logistic regression models):
agg_otu <- agg_otu_data %>% 
  mutate(key=otu) %>% 
  group_by(key)
taxa_info <- read.delim('data/process/vendors.taxonomy', header=T, sep='\t') %>% 
  select(-Size) %>% 
  mutate(key=OTU) %>% 
  select(-OTU)
agg_otu_data <- inner_join(agg_otu, taxa_info, by="key") %>%
  ungroup() %>% 
  mutate(key=str_to_upper(key)) %>% 
  mutate(taxa=gsub("(.*);.*","\\1",Taxonomy)) %>% 
  mutate(taxa=gsub("(.*)_.*","\\1",Taxonomy)) %>% 
  mutate(taxa=gsub("(.*);.*","\\1",Taxonomy)) %>% 
  mutate(taxa=gsub(".*;","",taxa)) %>% 
  mutate(taxa=gsub("(.*)_.*","\\1",taxa)) %>% 
  mutate(taxa=gsub('[0-9]+', '', taxa)) %>% 
  mutate(taxa=str_remove_all(taxa, "[(100)]")) %>% 
  unite(key, taxa, key, sep=" (") %>% 
  mutate(key = paste(key,")", sep="")) %>% 
  select(-otu, -Taxonomy) %>% 
  rename(otu=key) %>% 
  mutate(otu=paste0(gsub('TU0*', 'TU ', otu))) %>% 
  separate(otu, into = c("bactname", "OTUnumber"), sep = "\\ [(]", remove = FALSE) %>% #Add columns to separate bacteria name from OTU number to utilize ggtext so that only bacteria name is italicized
  mutate(otu_name = glue("*{bactname}* ({OTUnumber}")) #Markdown notation so that only bacteria name is italicized

#Hypothesis testing & plotting----

#Function to pull significant taxa (adjusted p value < 0.05) after statistical analysis
pull_significant_taxa <- function(dataframe, taxonomic_level){
  dataframe %>% 
    filter(p.value.adj <= 0.05) %>% 
    pull({{ taxonomic_level }}) #Embracing transforms taxonomic_level argument into a column name
}

#List of days with sequence data----
exp_days_seq <- unique(agg_otu_data %>% pull(day))
  
#Function for Kruskal_wallis test for differences across groups at different taxonomic levels with Benjamini-Hochburg correction----
set.seed(19881117) #Same seed used for mothur analysis

#Function to test at the otu level:
kruskal_wallis_otu <- function(timepoint){
  otu_stats <- agg_otu_data %>% 
    filter(day == timepoint) %>%
    select(vendor, otu, agg_rel_abund) %>% 
    group_by(otu) %>% 
    nest() %>% 
    mutate(model=map(data, ~kruskal.test(x=.x$agg_rel_abund, g=as.factor(.x$vendor)) %>% tidy())) %>% 
    mutate(median = map(data, get_rel_abund_median_vendor)) %>% 
    unnest(c(model, median)) %>% 
    ungroup() 
  #Adjust p-values for testing multiple OTUs
  otu_stats_adjust <- otu_stats %>% 
    select(otu, statistic, p.value, parameter, method, Schloss, Young, Jackson, `Charles River`, Taconic, Envigo) %>% 
    mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
    arrange(p.value.adj) %>% 
    write_tsv(path = paste0("data/process/otu_stats_day_", timepoint, ".tsv"))
}

# Perform kruskal wallis tests at the otu level for all days of the experiment that were sequenced----
for (d in exp_days_seq){
  kruskal_wallis_otu(d)
  #Make a list of significant otus across sources of mice for a specific day  
  stats <- read_tsv(file = paste0("data/process/otu_stats_day_", d, ".tsv"))
  name <- paste("sig_otu_day", d, sep = "") 
  assign(name, pull_significant_taxa(stats, otu))
}

#Shared significant genera across D-1 to D1----
shared_sig_otus_Dn1toD1 <- intersect_all(`sig_otu_day-1`, sig_otu_day0, sig_otu_day1)
# 12 OTUs: [1] "Parasutterella (OTU 26)" "Burkholderiales (OTU 34)" "Betaproteobacteria (OTU 58)" "Lactobacillus (OTU 49)"     
# "Parabacteroides (OTU 5)" "Lactobacillus (OTU 31)" "Bacteroides (OTU 2)" "Lachnospiraceae (OTU 130)"  
# "Proteus (OTU 16)" "Lactobacillus (OTU 6)" "Enterococcus (OTU 23)" "Ruminococcaceae (OTU 152)"   

#Function to plot a list of OTUs across sources of mice at a specific timepoint:
#Arguments: otus = list of otus to plot; timepoint = day of the experiment to plot
plot_otus_dx <- function(otus, timepoint){
  agg_otu_data %>% 
    filter(otu %in% otus) %>% 
    filter(day == timepoint) %>% 
    mutate(agg_rel_abund = agg_rel_abund + 1/10874) %>% # 10,874 is 2 times the subsampling parameter of 5437
    ggplot(aes(x= otu_name, y=agg_rel_abund, color=vendor))+
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_vendors,
                        labels=color_vendors)+
    geom_hline(yintercept=1/5437, color="gray")+
    stat_summary(fun = 'median', 
                 fun.max = function(x) quantile(x, 0.75), 
                 fun.min = function(x) quantile(x, 0.25),
                 position = position_dodge(width = 1)) +
    labs(title=NULL, 
         x=NULL,
         y="Relative abundance (%)")+
    scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100), limits = c(1/10900, 1))+
    coord_flip()+
    theme_classic()+
    theme(plot.title=element_text(hjust=0.5),
          legend.position = "bottom",
          axis.text.y = element_markdown(), #Have only the OTU names show up as italics
          text = element_text(size = 19)) # Change font size for entire plot
}

#Plots of the top 18-20 OTUs that varied across sources at each timepoint (only 18 OTUs varied significantly across sources post-clindamycin treatment)
#Baseline: top 20 OTUs that vary across sources
Dn1top20_otus <- plot_otus_dx(`sig_otu_day-1`[1:20], -1) +#Pick top 20 significant OTUs
  geom_vline(xintercept = c((1:20) - 0.5 ), color = "grey") + # Add gray lines to clearly separate OTUs
  theme(legend.position = "none") #remove legend
save_plot("results/figures/Dn1top20_otus.png", Dn1top20_otus, base_height = 9, base_width = 7)
#Post-clindamycin: top 18 OTUs that vary across sources
D0top18_otus <- plot_otus_dx(`sig_otu_day0`[1:18], 0) + #Pick top 20 significant OTUs
  geom_vline(xintercept = c((1:18) - 0.5 ), color = "grey") + # Add gray lines to clearly separate OTUs
  theme(legend.position = "none") #remove legend
save_plot("results/figures/D0top18_otus.png", D0top18_otus, base_height = 9, base_width = 7)
#Post-infection: top 20 OTUs that vary across sources
D1top20_otus <- plot_otus_dx(`sig_otu_day1`[1:20], 1)+ #Pick top 20 significant OTUs
  geom_vline(xintercept = c((1:18) - 0.5 ), color = "grey") + # Add gray lines to clearly separate OTUs
  theme(legend.position = "none") #remove legend
save_plot("results/figures/D1top20_otus.png", D1top20_otus, base_height = 9, base_width = 7)

#Overlap between the top 18-20 OTUs that varied across sources at each timepoints
shared_top_otus_Dn1toD1 <- intersect_all(`sig_otu_day-1`[1:20], sig_otu_day0[1:18], sig_otu_day1[1:20])
#5 OTUs: 1] "Parasutterella (OTU 26)", "Burkholderiales (OTU 34)", "Betaproteobacteria (OTU 58)","Lactobacillus (OTU 49)","Parabacteroides (OTU 5)"  

# Perform pairwise Wilcoxan rank sum tests for otus that were significantly different across sources of mice on a series of days----
pairwise_day_otu <- function(timepoint, sig_otu_dayX){
  otu_stats <- agg_otu_data %>% 
    filter(day == timepoint) %>%
    select(vendor, otu, agg_rel_abund) %>% 
    group_by(otu) %>% 
    nest() %>% 
    mutate(model=map(data, ~kruskal.test(x=.x$agg_rel_abund, g=as.factor(.x$vendor)) %>% tidy())) %>% 
    mutate(median = map(data, get_rel_abund_median_vendor)) %>% 
    unnest(c(model, median)) %>% 
    ungroup()
  pairwise_stats <- otu_stats %>% 
    filter(otu %in% sig_otu_dayX) %>% 
    group_by(otu) %>% 
    mutate(model=map(data, ~pairwise.wilcox.test(x=.x$agg_rel_abund, g=as.factor(.x$vendor), p.adjust.method="BH") %>% 
                       tidy() %>% 
                       mutate(compare=paste(group1, group2, sep="-")) %>% 
                       select(-group1, -group2) %>% 
                       pivot_wider(names_from=compare, values_from=p.value)
    )
    ) %>% 
    unnest(model) %>% 
    select(-data, -parameter, -statistic, -p.value) %>% #Get rid of p.value since it's the unadjusted version
    write_tsv(path = paste0("data/process/otu_stats_day_", timepoint, "_sig.tsv"))
  #Format pairwise stats to use with ggpubr package
  plot_format_stats <- pairwise_stats %>% 
    select(-method,-Schloss, -Young, -Jackson, -`Charles River`, -Taconic, -Envigo) %>% 
    group_split() %>% #Keeps a attr(,"ptype") to track prototype of the splits
    lapply(tidy_pairwise_otu) %>% 
    bind_rows()
  return(plot_format_stats)  
}
#Perform pairwise comparisons for day -1, 0, and 1. When mice are undergoing greatest perturbations (clindamycin administration, followed by C. difficile challenge)
otu_dayn1_stats <- pairwise_day_otu(-1, `sig_otu_day-1`)
otu_day0_stats <- pairwise_day_otu(0, sig_otu_day0)
otu_day1_stats <- pairwise_day_otu(1, sig_otu_day1)

#Combine tables of Kruskal-Wallis and Wilcoxan rank sum pairwise comparisons between sources of mice for day -1 to day 0 (3 timepoints used as input data for logistic regression models----
#Pull in day -1 tables:
w_sig_dn1 <- read_tsv(file="data/process/otu_stats_day_-1_sig.tsv") %>%
  select(-method, -Schloss, -Young, -Jackson, -`Charles River`, -Taconic, -Envigo) #Remove duplicate columns
#Combine Kruskal-Wallis and Wilcoxan rank sum pairwise for each day
kw_w_sig_dn1 <- read_tsv(file="data/process/otu_stats_day_-1.tsv") %>%  
  select(-parameter) %>% #Don't need this column
  inner_join(w_sig_dn1, by = c("otu")) %>% 
  mutate(day = -1) %>% #Add experiment day for statistics
  arrange(p.value.adj) #Arrange by adjusted p.values

#Pull in day 0 tables:
w_sig_d0 <- read_tsv(file="data/process/otu_stats_day_0_sig.tsv") %>%
  select(-method, -Schloss, -Young, -Jackson, -`Charles River`, -Taconic, -Envigo) #Remove duplicate columns
#Combine Kruskal-Wallis and Wilcoxan rank sum pairwise for each day
kw_w_sig_d0 <- read_tsv(file="data/process/otu_stats_day_0.tsv") %>%  
  select(-parameter) %>% #Don't need this column
  inner_join(w_sig_d0, by = c("otu")) %>% 
  mutate(day = 0) %>% #Add experiment day for statistics
  arrange(p.value.adj) #Arrange by adjusted p.values

#Pull in day 1 tables:
w_sig_d1 <- read_tsv(file="data/process/otu_stats_day_1_sig.tsv") %>%
  select(-method, -Schloss, -Young, -Jackson, -`Charles River`, -Taconic, -Envigo) #Remove duplicate columns
#Combine Kruskal-Wallis and Wilcoxan rank sum pairwise for each day
kw_w_sig_d1 <- read_tsv(file="data/process/otu_stats_day_1.tsv") %>%  
  select(-parameter) %>% #Don't need this column
  inner_join(w_sig_d1, by = c("otu")) %>% 
  mutate(day = 1) %>% #Add experiment day for statistics
  arrange(p.value.adj) #Arrange by adjusted p.values

#Combine OTU statistics for differences across sources of mice for days -1, 0 and 1:
otu_stats_dn1to1_combined <- rbind(kw_w_sig_dn1, kw_w_sig_d0, kw_w_sig_d1) %>% 
  write_tsv(path = paste0("data/process/otu_stats_dn1to1_combined.tsv"))  #Save combined dataframe as a .tsv

#Read in just Kruskal Wallis results for the rest of the days:
kw_sig_dn1 <- read_tsv(file="data/process/otu_stats_day_-1.tsv") %>%  
  select(-parameter) %>% #Don't need this column
  mutate(day = -1) %>% #Add experiment day for statistics
  arrange(p.value.adj) #Arrange by adjusted p.values
kw_sig_d0 <- read_tsv(file="data/process/otu_stats_day_0.tsv") %>%  
  select(-parameter) %>% #Don't need this column
  mutate(day = 0) %>% #Add experiment day for statistics
  arrange(p.value.adj) #Arrange by adjusted p.values
kw_sig_d1 <- read_tsv(file="data/process/otu_stats_day_1.tsv") %>%  
  select(-parameter) %>% #Don't need this column
  mutate(day = 1) %>% #Add experiment day for statistics
  arrange(p.value.adj) #Arrange by adjusted p.values
kw_sig_d2 <- read_tsv(file="data/process/otu_stats_day_2.tsv") %>%  
  select(-parameter) %>% #Don't need this column
  mutate(day = 2) %>% #Add experiment day for statistics
  arrange(p.value.adj) #Arrange by adjusted p.values
kw_sig_d3 <- read_tsv(file="data/process/otu_stats_day_3.tsv") %>%  
  select(-parameter) %>% #Don't need this column
  mutate(day = 3) %>% #Add experiment day for statistics
  arrange(p.value.adj) #Arrange by adjusted p.values
kw_sig_d4 <- read_tsv(file="data/process/otu_stats_day_4.tsv") %>%  
  select(-parameter) %>% #Don't need this column
  mutate(day = 4) %>% #Add experiment day for statistics
  arrange(p.value.adj) #Arrange by adjusted p.values
kw_sig_d5 <- read_tsv(file="data/process/otu_stats_day_5.tsv") %>%  
  select(-parameter) %>% #Don't need this column
  mutate(day = 5) %>% #Add experiment day for statistics
  arrange(p.value.adj) #Arrange by adjusted p.values
kw_sig_d6 <- read_tsv(file="data/process/otu_stats_day_6.tsv") %>%  
  select(-parameter) %>% #Don't need this column
  mutate(day = 6) %>% #Add experiment day for statistics
  arrange(p.value.adj) #Arrange by adjusted p.values
kw_sig_d7 <- read_tsv(file="data/process/otu_stats_day_7.tsv") %>%  
  select(-parameter) %>% #Don't need this column
  mutate(day = 7) %>% #Add experiment day for statistics
  arrange(p.value.adj) #Arrange by adjusted p.values
kw_sig_d8 <- read_tsv(file="data/process/otu_stats_day_8.tsv") %>%  
  select(-parameter) %>% #Don't need this column
  mutate(day = 8) %>% #Add experiment day for statistics
  arrange(p.value.adj) #Arrange by adjusted p.values
kw_sig_d9 <- read_tsv(file="data/process/otu_stats_day_9.tsv") %>%  
  select(-parameter) %>% #Don't need this column
  mutate(day = 9) %>% #Add experiment day for statistics
  arrange(p.value.adj) #Arrange by adjusted p.values
#Combine Kruskal-Wallis results at the OTU level for all days with sequence data:----
otu_stats_dn1to9_combined <- rbind(kw_sig_dn1, kw_sig_d0, kw_sig_d1, kw_sig_d2, kw_sig_d3, kw_sig_d4,
                                      kw_sig_d5, kw_sig_d6, kw_sig_d7, kw_sig_d8, kw_sig_d9) %>% 
  filter(p.value.adj <= 0.05) %>% #Select only significant rows
  write_tsv(path = paste0("data/process/otu_kw_stats_dn1to9.tsv")) #Save combined dataframe as a .tsv

#Wilcoxan Signed rank test for relative abundance differences after clindamycin treatment (for all mice with paired data for day -1 versus day 0) at different taxonomic levels with Benjamini-Hochburg correction----
#Pull mice ids that have sequence data for both day -1 and day 0:
mice_seq_dn1_0_pairs <- agg_otu_data %>% 
  filter(otu == "Bacteroides (OTU 2)") %>% #Random pick just to figure out what mice have sequence data
  filter(day == -1 | day == 0) %>% 
  filter(duplicated(mouse_id)) %>% #Pull mouse ids with sequence data for both day -1 and day 0
  pull(mouse_id) 

#Dataframe for statistical test at the OTU level
paired_otu <- agg_otu_data %>% 
  filter(mouse_id %in% mice_seq_dn1_0_pairs) %>% #Only select pairs with data for day -1 & day 0
  filter(day == -1 | day == 0) %>% #Experiment days that represent initial community and community post clindamycin treatment
  mutate(day = as.factor(day)) %>% 
  select(day, otu, agg_rel_abund)

#Wilcoxon signed rank test for all day -1, day 0 pairs at the OTU level:
o_dn1to0_pairs <- paired_otu %>% 
  group_by(otu) %>% 
    nest() %>% 
    mutate(model=map(data, ~wilcox.test(.x$agg_rel_abund ~ .x$day, paired = TRUE) %>% tidy())) %>% 
    mutate(median = map(data, get_rel_abund_median_day)) %>% 
    unnest(c(model, median)) %>% 
    ungroup() 
#Adjust p-values for testing multiple OTUs
o_dn1to0_pairs_stats_adjust <- o_dn1to0_pairs %>% 
    select(otu, statistic, p.value, method, alternative, `-1`, `0`) %>% 
    mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
    arrange(p.value.adj) %>% 
  write_tsv(path = "data/process/otu_stats_dn1to0.tsv") 

#Make a list of significant OTUs impacted by clindamycin treatment----  
sig_otu_pairs <- pull_significant_taxa(o_dn1to0_pairs_stats_adjust, otu)
# 153 OTUs
sig_otu_pairs_top10 <- sig_otu_pairs[1:10]

#OTUs that overlap between Fig. 3D and E----
Fig3D_vs_E_overlap <- intersect_all(sig_otu_day0[1:18], sig_otu_pairs[1:10])
#"Enterobacteriaceae (OTU 1)"
Clind_vs_source_D0 <- intersect_all(sig_otu_day0, sig_otu_pairs)
#4 OTUs: "Enterobacteriaceae (OTU 1)","Enterococcus (OTU 23)","Lactobacillus (OTU 6)","Lachnospiraceae (OTU 130)"

#Plot of the otus with significantly different relative abundances post clindamycin treatment across all sources of mice:----
#Facet plot by day with the following label names:
facet_labels <- c("Baseline", "Clindamycin")
names(facet_labels) <- c("-1", "0")
clind_impacted_otus_plot_dn1_0 <- agg_otu_data %>% 
  filter(mouse_id %in% mice_seq_dn1_0_pairs) %>% #Only select pairs with data for day -1 & day 0
  filter(otu %in% sig_otu_pairs_top10) %>% 
  filter(day %in% c(-1, 0)) %>% #Select baseline and post-clindamycin timepoints
  mutate(otu=factor(otu, levels=sig_otu_pairs_top10)) %>% 
  mutate(agg_rel_abund = agg_rel_abund + 1/10874) %>% # 10,874 is 2 times the subsampling parameter of 5437
  ggplot(aes(x= otu_name, y=agg_rel_abund, color=vendor))+
  scale_colour_manual(name=NULL,
                      values=color_scheme,
                      breaks=color_vendors,
                      labels=color_vendors)+
  geom_hline(yintercept=1/5437, color="gray")+
    stat_summary(fun = 'median', 
                 fun.max = function(x) quantile(x, 0.75), 
                 fun.min = function(x) quantile(x, 0.25),
                 position = position_dodge(width = 1)) +  
  labs(title=NULL, 
       x=NULL,
       y="Relative abundance (%)")+
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100), limits = c(1/10900, 1))+
  coord_flip()+
  theme_classic()+
  geom_vline(xintercept = c((1:10) - 0.5 ), color = "grey") + # Add gray lines to clearly separate OTUs
  facet_wrap( ~ day, labeller = labeller(day = facet_labels), scales = "fixed")+
  theme(plot.title=element_text(hjust=0.5),
        text = element_text(size = 19),# Change font size for entire plot
        axis.text.y = element_markdown(), #Have only the OTU names show up as italics
        strip.background = element_blank(),
        legend.position = "none") 
save_plot(filename = paste0("results/figures/clind_impacted_otus_plot.png"), clind_impacted_otus_plot_dn1_0, base_height = 9, base_width = 7)

#Comparison of Figure 4 (varied across colony source) and 5 (altered by clindamycin treatment) taxa----
Fig4_v_5_otus <- intersect_all(`shared_sig_otus_Dn1toD1`, `sig_otu_pairs_top10`)
#O OTUs overlap
Fig4_v_all_clind <- intersect_all(`shared_sig_otus_Dn1toD1`, `sig_otu_pairs`)
# 3 OTUs Lachnospiraceae (OTU 130), Lactobacillus (OTU 6), Enterococcus (OTU 23)

#Function to plot OTUs of interest that overlap with top 20 OTUs in 3 logistic regression models----
#Customize the x_annotation, y_position, and label argument values for each OTU prior to running the function
otu_over_time <- function(otu_plot, x_annotation, y_position, label){
  specify_otu_name <- agg_otu_data %>% 
    filter(otu == otu_plot) %>% 
    pull(otu_name)
  otu_median <- agg_otu_data %>% 
    filter(otu == otu_plot) %>% 
    group_by(vendor, day) %>% 
    summarize(median=(median(agg_rel_abund + 1/10874))) %>% 
    ungroup
  otu_mice <- agg_otu_data %>% 
    filter(otu == otu_plot) %>% 
    mutate(agg_rel_abund = agg_rel_abund + 1/10874) %>%
    select(vendor, day, agg_rel_abund, otu, experiment, clearance_status_d7)
  otu_time <- ggplot(NULL)+
    geom_point(otu_mice, mapping = aes(x=day, y=agg_rel_abund, color=vendor, shape=clearance_status_d7), size  = 1.5, position = position_dodge(width = 0.6))+
    geom_line(otu_median, mapping = aes(x=day, y=median, color=vendor), size = 1, show.legend = FALSE)+
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_vendors,
                        labels=color_vendors)+
    scale_shape_manual(name="Cleared by Day 7",
                       values=c(4, 19, 21),
                       breaks=c("colonized", "not_detectable", "no_data"),
                       labels=c("no", "yes", "no data"), 
                       drop=FALSE, na.translate = TRUE, na.value = 1)+
    geom_hline(yintercept=1/5437, color="gray")+
    labs(title=specify_otu_name,
         x="Day",
         y="Relative abundance (%)") +
    scale_x_continuous(breaks = c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
                       limits = c(-1.5, 9.5)) +
    scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
    annotate("text", y = y_position, x = x_annotation, label = label, size =8)+
    theme_classic()+
    theme(plot.title=element_markdown(hjust = 0.5),
          panel.grid.minor.x = element_line(size = 0.4, color = "grey"),  # Add gray lines to clearly separate symbols by days)
          text = element_text(size = 18)) # Change font size for entire plot
}

#OTUs that differ across mouse sources at multiple timepoints and are important features in at least 2/3 logistic regresssion models----
#Create column with statistical significance symbols for all OTUs:
plot_otu_stats_dn1to9 <- otu_stats_dn1to9_combined %>% 
  mutate(p.value.adj=round(p.value.adj, digits = 4)) %>% 
  filter(p.value.adj <= 0.05) %>%
  mutate(p.signif = case_when(
    p.value.adj > 0.05 ~ "NS",
    p.value.adj <= 0.05 ~ "*"
  )) 
#Set up statistical annotation arguments for Bacteroides (OTU 2):
x_OTU2 <- plot_otu_stats_dn1to9 %>% 
  filter(otu == "Bacteroides (OTU 2)") %>% pull(day)
y_OTU2 <- 1.5
label_OTU2 <- plot_otu_stats_dn1to9 %>% 
  filter(otu == "Bacteroides (OTU 2)") %>% pull(p.signif)
otu2 <- otu_over_time(otu_plot = "Bacteroides (OTU 2)", x_annotation = x_OTU2, y_position = y_OTU2, label = label_OTU2)+
  theme(legend.position = "none")
save_plot(filename = paste0("results/figures/Bacteroides (OTU 2)_time.png"), otu2, base_aspect_ratio = 2.5)

#OTUs that vary across sources, are impacted by clindamycin and were important in at least 2/3 logistic regression models:
#Set up statistical annotation arguments for Enterococcus (OTU 23):
x_OTU23 <- plot_otu_stats_dn1to9 %>% 
  filter(otu == "Enterococcus (OTU 23)") %>% pull(day)
y_OTU23 <- .4
label_OTU23 <- plot_otu_stats_dn1to9 %>% 
  filter(otu == "Enterococcus (OTU 23)") %>% pull(p.signif)
otu23 <- otu_over_time("Enterococcus (OTU 23)", x_annotation = x_OTU23, y_position = y_OTU23, label = label_OTU23)+
  theme(legend.position = "none")
save_plot(filename = paste0("results/figures/Enterococcus (OTU 23)_time.png"), otu23, base_aspect_ratio = 2.5)

#Set up statistical annotation arguments for Enterobacteriaceae (OTU 1):
x_OTU1 <- plot_otu_stats_dn1to9 %>% 
  filter(otu == "Enterobacteriaceae (OTU 1)") %>% pull(day)
y_OTU1 <- 1.5
label_OTU1 <- plot_otu_stats_dn1to9 %>% 
  filter(otu == "Enterobacteriaceae (OTU 1)") %>% pull(p.signif)
otu1 <- otu_over_time("Enterobacteriaceae (OTU 1)", x_annotation = x_OTU1, y_position = y_OTU1, label = label_OTU1)+
  theme(legend.position = "none")
save_plot(filename = paste0("results/figures/Enterobacteriaceae (OTU 1)_time.png"), otu1, base_aspect_ratio = 2.5)

#OTUs impacted by clindamycin treatment and are important features in at least 2/3 logistic regresssion models---- 
#Set up statistical annotation arguments for Porphyromonadaceae (OTU 7):
x_OTU7 <- plot_otu_stats_dn1to9 %>% 
  filter(otu == "Porphyromonadaceae (OTU 7)") %>% pull(day)
y_OTU7 <- .4
label_OTU7 <- plot_otu_stats_dn1to9 %>% 
  filter(otu == "Porphyromonadaceae (OTU 7)") %>% pull(p.signif)
otu7 <- otu_over_time("Porphyromonadaceae (OTU 7)", x_annotation = x_OTU7, y_position = y_OTU7, label = label_OTU7)+
  theme(legend.position = "bottom", 
        legend.text = element_text(size = 16))
legend <- get_legend(otu7)
#Alternative function to get legend:
legend <- cowplot::get_legend(otu7)
overlap_OTUs_legend <- as_ggplot(legend)#+theme(element_line(size = 16)) #+theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
save_plot(filename = paste0("results/figures/overlap_OTUs_legend.png"), overlap_OTUs_legend, base_height = .9, base_width = 9.5)
 otu7 <- otu7 + theme(legend.position = "none") #remove legend
save_plot(filename = paste0("results/figures/Porphyromonadaceae (OTU 7)_time.png"), otu7, base_aspect_ratio = 2.5)

#Examine which taxa have different relative abundances between experiments for Schloss, Young, and Envigo mice at baseline----
#Function to test for differences in relative abundances at the otu level according to experiment within specific sources of mice:
baseline_exp_o <- function(vendor_name){
  otu_stats <- agg_otu_data %>% 
    filter(day == -1) %>% 
    filter(vendor == vendor_name) %>% 
    select(experiment, otu, agg_rel_abund) %>% 
    group_by(otu) %>% 
    nest() %>% 
    mutate(model=map(data, ~wilcox.test(.x$agg_rel_abund ~ .x$experiment, paired = FALSE) %>% tidy())) %>% 
    mutate(median = map(data, get_rel_abund_median_experiment)) %>% 
    unnest(c(model, median)) %>% 
    ungroup() 
  #Adjust p-values for testing multiple OTUs
  otu_stats_adjust <- otu_stats %>% 
    select(otu, statistic, p.value, method, alternative, `1`, `2`) %>% 
    mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
    arrange(p.value.adj) %>% 
    write_tsv(path = paste0("data/process/otu_stats_experiment_", vendor_name, ".tsv"))
}

#Test only Schloss, Young, and Envigo mice
test_sources <- c("Schloss", "Young", "Envigo")
for (v in test_sources){
  baseline_exp_o(v)
  #Make a list of significant OTUs that differ according to experiment for a specific source of mice 
  stats <- read_tsv(file = paste0("data/process/otu_stats_experiment_", v, ".tsv"))
  name <- paste("sig_otu_experiment_", v, sep = "") 
  assign(name, pull_significant_taxa(stats, otu))
}

sig_otu_experiment_Schloss #No significant OTUs
sig_otu_experiment_Young #No significant OTUs
sig_otu_experiment_Envigo #No significant OTUs

#Since no OTUs were significant after multiple hypothesis correction. Check if any OTUs had p.value < 0.05
#Function to pull different taxa (p value < 0.05) before BH adjustment
pull_different_taxa <- function(dataframe, taxonomic_level){
  dataframe %>% 
    filter(p.value <= 0.05) %>% 
    pull({{ taxonomic_level }}) #Embracing to transform of taxonomic_level argument into a column name
}
test_sources <- c("Schloss", "Young", "Envigo")
for (v in test_sources){
  #Make a list of OTUs that differ according to experiment for a specific source of mice and had p.value < 0.05 before BH adjustment
  stats <- read_tsv(file = paste0("data/process/otu_stats_experiment_", v, ".tsv"))
  name <- paste("diff_otu_experiment_", v, sep = "") 
  assign(name, pull_different_taxa(stats, otu))
}

diff_otu_experiment_Schloss #0
diff_otu_experiment_Young #58
diff_otu_experiment_Envigo #3
  #Porphyromonadaceae OTU 52 and 139, Deltaproteobacteria OTU 174

#Check if any overlap with taxa that varied between sources at baseline
Young_exp_baseline <- intersect_all(diff_otu_experiment_Young, `sig_otu_day-1`)
#21 overlap
Young_exp_baseline <- intersect_all(diff_otu_experiment_Young, `sig_otu_day-1`[1:20])
#O overlap with top 20
Envigo_exp_baseline <- intersect_all(diff_otu_experiment_Envigo, `sig_otu_day-1`)
#Deltaproteobacteria (OTU 174) 

#Examine top 20 taxa from logistic regression model that performed the best: Communities after clindamycin treatment at the OTU level----
interp_otus <- read_tsv("data/process/combined_top20_otus_all_models.tsv")

#List of top 20 taxa from day 0 OTU logistic regression model
interp_otus_d0 <- interp_otus  %>% filter(model_input_day == 0) %>% pull(OTU)

#Create a data frame of just the top 10 most important OTUs for the Day 0 logistic regression model
#Create a color column based on correlation coefficient sign.
interp_otus_d0_top_10 <- interp_otus  %>% filter(model_input_day == 0) %>% 
  arrange(median_rank) %>% 
  slice(1:10) %>% #Pick top 10 OTUs based on median_rank
  #Assign color to OTU based on correlation coefficient (- correlates with colonization, red; + correlates with clearance, dark blue)
  mutate(color = case_when(median_feature_weight < 0 ~ "firebrick",
                           median_feature_weight > 0 ~ "navyblue"))
  
#Plot the top 10 OTUs that were important to Day 0 model and 
#using facet_wrap, highlight the OTUs that correlate with colonization
top10_d0_otu_model_taxa <- agg_otu_data %>% 
  filter(day == 0) %>% 
  filter(otu %in% (interp_otus_d0_top_10 %>% pull(OTU))) %>%
  mutate(agg_rel_abund = agg_rel_abund + 1/10874) %>% 
  group_by(vendor, otu) %>% 
  left_join(interp_otus_d0_top_10, by = c("otu" = "OTU")) %>%  #Join to interp_otus_d0_top_10 to preserve coeffic.
  #Modify ggtext formatted OTU label to incorp. color based on correlation coefficient sign in d0_model
  mutate(otu_name = glue("<span style='color:{color}'>*{bactname}* ({OTUnumber}")) %>%  #Markdown notation so that only bacteria name is italicized
  mutate(median=(median(agg_rel_abund))) %>% #create a column of median values for each group
  ungroup() %>% 
  arrange(median_rank) %>% 
  ggplot(aes(x=vendor, y =agg_rel_abund, colour= vendor))+
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  geom_errorbar(aes(ymax = median, ymin = median), color = "gray50", size = 1)+ #Add lines to indicate the median for each group to the plot
  geom_jitter(aes(shape = clearance_status_d7), size=2, show.legend = TRUE) +
  scale_colour_manual(name=NULL,
                      values=color_scheme,
                      breaks=color_vendors,
                      labels=color_vendors)+
  scale_shape_manual(name="Cleared by Day 7",
                     values=c(4, 19, 21),
                     breaks=c("colonized", "not_detectable", "no_data"),
                     labels=c("no", "yes", "no data"), 
                     drop=FALSE, na.translate = TRUE, na.value = 1)+
  geom_hline(yintercept=1/5437, color="gray")+
  facet_wrap(~ otu_name, nrow=2)+
  labs(title=NULL,
       x=NULL,
       y="Relative abundance (%)") +
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
  theme_classic()+
  theme(text = element_text(size = 14),
        strip.text = element_markdown(hjust = 0.5, size = 6.8),
        axis.text.x = element_blank(),
        legend.position = "bottom")
save_plot(filename = paste0("results/figures/d0_top10_otus.png"), top10_d0_otu_model_taxa, base_height = 5, base_width = 8)


