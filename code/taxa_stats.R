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
  mutate(otu=paste0(gsub('TU0*', 'TU ', otu))) 

# Relative abundance data at the family level:
agg_family_data <- agg_taxonomic_data(family)

#Hypothesis testing & plotting----

#Function to pull significant taxa (adjusted p value < 0.05) after statistical analysis
pull_significant_taxa <- function(dataframe, taxonomic_level){
  taxonomic_level <- enquo(taxonomic_level) #Part of transformation of taxonomic_level argument into a column name
  dataframe %>% 
    filter(p.value.adj <= 0.05) %>% 
    pull(!!taxonomic_level) #!!Completes the transformation of taxonomic_level argument into a column name
}

#List of days with sequence data----
exp_days_seq <- unique(agg_otu_data %>% pull(day))
  
#Function for Kruskal_wallis test for differences across groups at different taxonomic levels with Benjamini-Hochburg correction----
set.seed(19881117) #Same seed used for mothur analysis
#Function to test at the family level:
kruskal_wallis_f <- function(timepoint){
  family_stats <- agg_family_data %>% 
  filter(day == timepoint) %>%
  select(vendor, family, agg_rel_abund) %>% 
  group_by(family) %>% 
  nest() %>% 
  mutate(model=map(data, ~kruskal.test(x=.x$agg_rel_abund, g=as.factor(.x$vendor)) %>% tidy())) %>% 
  mutate(median = map(data, get_rel_abund_median_vendor)) %>% 
  unnest(c(model, median)) %>% 
  ungroup() 
  #Adjust p-values for testing multiple families
  family_stats_adjust <- family_stats %>% 
    select(family, statistic, p.value, parameter, method, Schloss, Young, Jackson, `Charles River`, Taconic, Envigo) %>% 
    mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
    arrange(p.value.adj) %>% 
    write_tsv(path = paste0("data/process/family_stats_day_", timepoint, ".tsv"))
}

# Perform kruskal wallis tests at the family level for all days of the experiment that were sequenced----
for (d in exp_days_seq){
  kruskal_wallis_f(d)
  #Make a list of significant families across sources of mice for a specific day  
  stats <- read_tsv(file = paste0("data/process/family_stats_day_", d, ".tsv"))
  name <- paste("sig_family_day", d, sep = "") 
  assign(name, pull_significant_taxa(stats, family))
}

#Shared significant familes across D-1 to D1----
shared_sig_families_Dn1toD1 <- intersect_all(`sig_family_day-1`, sig_family_day0, sig_family_day1)
# 8 families: Betaproteobacteria Unclassified, Burkholderiales Unclassified, Sutterellaceae, Deferribacteraceae, Bacteroidaceae

#Function to plot a list of families across sources of mice at a specific timepoint:
#Arguments: families = list of families to plot; timepoint = day of the experiment to plot
plot_families_dx <- function(families, timepoint){
  agg_family_data %>% 
  filter(family %in% families) %>% 
  filter(day == timepoint) %>% 
  mutate(family=factor(family, levels=families)) %>% 
  mutate(agg_rel_abund = agg_rel_abund + 1/10874) %>% # 10,874 is 2 times the subsampling parameter of 5437
#  group_by(vendor, family) %>%   
#  mutate(median = median(agg_rel_abund)) %>% #create a column of median_cfu
#  ungroup() %>%      
  ggplot(aes(x = family, y=agg_rel_abund, color=vendor))+
  scale_colour_manual(name=NULL,
                      values=color_scheme,
                      breaks=color_vendors,
                      labels=color_vendors)+
  geom_hline(yintercept=1/5437, color="gray")+
  geom_boxplot(outlier.shape = NA, size = 1.2)+  
#  geom_errorbar(aes(ymax = median, ymin = median))+ #Add lines to indicate the median for each group to the plot. Median calculated before y axis transformation
#  geom_jitter(aes(shape = experiment), size=2, alpha=0.6, position=position_dodge(1.6)) +
#  scale_shape_manual(name=NULL,
#                     values=shape_scheme,
#                     breaks=shape_experiment,
#                     labels=shape_experiment) +
  labs(title=NULL, 
       x=NULL,
       y="Relative abundance (%)")+
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100), limits = c(1/10900, 1))+
  coord_flip()+
  theme_classic()+
  theme(plot.title=element_text(hjust=0.5),
        axis.text.y = element_text(face = "italic"), #Have the families show up as italics
        legend.position = "bottom",
        text = element_text(size = 16)) # Change font size for entire plot
}

#Plots of the relative abundances of families that significantly varied across sources of mice from day -1 to day 1----
Dn1toD1_families_dn1 <- plot_families_dx(shared_sig_families_Dn1toD1, -1) +
  ggtitle("Baseline")+ #Title plot
  theme(plot.title = element_text(hjust = 0.5)) #Center plot title
save_plot("results/figures/Dn1toD1_families_dn1.png", Dn1toD1_families_dn1, base_height = 6, base_width = 8)
Dn1toD1_families_d0 <- plot_families_dx(shared_sig_families_Dn1toD1, 0) +
  ggtitle("Clindamycin")+ #Title plot
  theme(plot.title = element_text(hjust = 0.5)) #Center plot title
save_plot("results/figures/Dn1toD1_families_d0.png", Dn1toD1_families_d0, base_height = 6, base_width = 8)
Dn1toD1_families_d1 <- plot_families_dx(shared_sig_families_Dn1toD1, 1) +
  ggtitle("Post-infection")+ #Title plot
  theme(plot.title = element_text(hjust = 0.5)) #Center plot title
save_plot("results/figures/Dn1toD1_families_d1.png", Dn1toD1_families_d1, base_height = 6, base_width = 8)

# Perform pairwise Wilcoxan rank sum tests for families that were significantly different across sources of mice on a series of days----
pairwise_day_family <- function(timepoint, sig_family_dayX){
  family_stats <- agg_family_data %>% 
    filter(day == timepoint) %>%
    select(vendor, family, agg_rel_abund) %>% 
    group_by(family) %>% 
    nest() %>% 
    mutate(model=map(data, ~kruskal.test(x=.x$agg_rel_abund, g=as.factor(.x$vendor)) %>% tidy())) %>% 
    mutate(median = map(data, get_rel_abund_median_vendor)) %>% 
    unnest(c(model, median)) %>% 
    ungroup()
  pairwise_stats <- family_stats %>% 
    filter(family %in% sig_family_dayX) %>% 
    group_by(family) %>% 
    mutate(model=map(data, ~pairwise.wilcox.test(x=.x$agg_rel_abund, g=as.factor(.x$vendor), p.adjust.method="BH") %>% 
                       tidy() %>% 
                       mutate(compare=paste(group1, group2, sep="-")) %>% 
                       select(-group1, -group2) %>% 
                       pivot_wider(names_from=compare, values_from=p.value)
    )
    ) %>% 
    unnest(model) %>% 
    select(-data, -parameter, -statistic, -p.value) %>% #Get rid of p.value since it's the unadjusted version
    write_tsv(path = paste0("data/process/family_stats_day_", timepoint, "_sig.tsv"))
  #Format pairwise stats to use with ggpubr package
  plot_format_stats <- pairwise_stats %>% 
    select(-method,-Schloss, -Young, -Jackson, -`Charles River`, -Taconic, -Envigo) %>% 
    group_split() %>% #Keeps a attr(,"ptype") to track prototype of the splits
    lapply(tidy_pairwise_family) %>% 
    bind_rows()
  return(plot_format_stats)  
}

#Perform pairwise comparisons for day -1, 0, and 1. When mice are undergoing greatest perturbations (clindamycin administration, followed by C. difficile challenge)
family_dayn1_stats <- pairwise_day_family(-1, `sig_family_day-1`)
family_day0_stats <- pairwise_day_family(0, sig_family_day0)
family_day1_stats <- pairwise_day_family(1, sig_family_day1)

#Combine tables of Kruskal-Wallis and Wilcoxan rank sum pairwise comparisons between sources of mice for day -1 to day 0 (3 timepoints used as input data for logistic regression models----
#Pull in day -1 tables:
w_sig_dn1 <- read_tsv(file="data/process/family_stats_day_-1_sig.tsv") %>%
  select(-method, -Schloss, -Young, -Jackson, -`Charles River`, -Taconic, -Envigo) #Remove duplicate columns
#Combine Kruskal-Wallis and Wilcoxan rank sum pairwise for each day
kw_w_sig_dn1 <- read_tsv(file="data/process/family_stats_day_-1.tsv") %>%  
  select(-parameter) %>% #Don't need this column
  inner_join(w_sig_dn1, by = c("family")) %>% 
  mutate(day = -1) %>% #Add experiment day for statistics
  arrange(p.value.adj) #Arrange by adjusted p.values

#Pull in day 0 tables:
w_sig_d0 <- read_tsv(file="data/process/family_stats_day_0_sig.tsv") %>%
  select(-method, -Schloss, -Young, -Jackson, -`Charles River`, -Taconic, -Envigo) #Remove duplicate columns
#Combine Kruskal-Wallis and Wilcoxan rank sum pairwise for each day
kw_w_sig_d0 <- read_tsv(file="data/process/family_stats_day_0.tsv") %>%  
  select(-parameter) %>% #Don't need this column
  inner_join(w_sig_d0, by = c("family")) %>% 
  mutate(day = 0) %>% #Add experiment day for statistics
  arrange(p.value.adj) #Arrange by adjusted p.values

#Pull in day 1 tables:
w_sig_d1 <- read_tsv(file="data/process/family_stats_day_1_sig.tsv") %>%
  select(-method, -Schloss, -Young, -Jackson, -`Charles River`, -Taconic, -Envigo) #Remove duplicate columns
#Combine Kruskal-Wallis and Wilcoxan rank sum pairwise for each day
kw_w_sig_d1 <- read_tsv(file="data/process/family_stats_day_1.tsv") %>%  
  select(-parameter) %>% #Don't need this column
  inner_join(w_sig_d1, by = c("family")) %>% 
  mutate(day = 1) %>% #Add experiment day for statistics
  arrange(p.value.adj) #Arrange by adjusted p.values

#Combine family statistics for differences across souces of mice for days -1, 0 and 1:----
family_stats_dn1to1_combined <- rbind(kw_w_sig_dn1, kw_w_sig_d0, kw_w_sig_d1) %>% 
  write_tsv(path = paste0("data/process/family_stats_dn1to1_combined.tsv")) %>% #Save combined dataframe as a .tsv
  #Also write results to supplemental table excel file
  write_xlsx("submission/table_S9_family_stats_dn1to1.xlsx", format_headers = FALSE)

#Read in just Kruskal Wallis results for the rest of the days:
kw_sig_dn1 <- read_tsv(file="data/process/family_stats_day_-1.tsv") %>%  
  select(-parameter) %>% #Don't need this column
  mutate(day = -1) %>% #Add experiment day for statistics
  arrange(p.value.adj) #Arrange by adjusted p.values
kw_sig_d0 <- read_tsv(file="data/process/family_stats_day_0.tsv") %>%  
  select(-parameter) %>% #Don't need this column
  mutate(day = 0) %>% #Add experiment day for statistics
  arrange(p.value.adj) #Arrange by adjusted p.values
kw_sig_d1 <- read_tsv(file="data/process/family_stats_day_1.tsv") %>%  
  select(-parameter) %>% #Don't need this column
  mutate(day = 1) %>% #Add experiment day for statistics
  arrange(p.value.adj) #Arrange by adjusted p.values
kw_sig_d2 <- read_tsv(file="data/process/family_stats_day_2.tsv") %>%  
  select(-parameter) %>% #Don't need this column
  mutate(day = 2) %>% #Add experiment day for statistics
  arrange(p.value.adj) #Arrange by adjusted p.values
kw_sig_d3 <- read_tsv(file="data/process/family_stats_day_3.tsv") %>%  
  select(-parameter) %>% #Don't need this column
  mutate(day = 3) %>% #Add experiment day for statistics
  arrange(p.value.adj) #Arrange by adjusted p.values
kw_sig_d4 <- read_tsv(file="data/process/family_stats_day_4.tsv") %>%  
  select(-parameter) %>% #Don't need this column
  mutate(day = 4) %>% #Add experiment day for statistics
  arrange(p.value.adj) #Arrange by adjusted p.values
kw_sig_d5 <- read_tsv(file="data/process/family_stats_day_5.tsv") %>%  
  select(-parameter) %>% #Don't need this column
  mutate(day = 5) %>% #Add experiment day for statistics
  arrange(p.value.adj) #Arrange by adjusted p.values
kw_sig_d6 <- read_tsv(file="data/process/family_stats_day_6.tsv") %>%  
  select(-parameter) %>% #Don't need this column
  mutate(day = 6) %>% #Add experiment day for statistics
  arrange(p.value.adj) #Arrange by adjusted p.values
kw_sig_d7 <- read_tsv(file="data/process/family_stats_day_7.tsv") %>%  
  select(-parameter) %>% #Don't need this column
  mutate(day = 7) %>% #Add experiment day for statistics
  arrange(p.value.adj) #Arrange by adjusted p.values
kw_sig_d8 <- read_tsv(file="data/process/family_stats_day_8.tsv") %>%  
  select(-parameter) %>% #Don't need this column
  mutate(day = 8) %>% #Add experiment day for statistics
  arrange(p.value.adj) #Arrange by adjusted p.values
kw_sig_d9 <- read_tsv(file="data/process/family_stats_day_9.tsv") %>%  
  select(-parameter) %>% #Don't need this column
  mutate(day = 9) %>% #Add experiment day for statistics
  arrange(p.value.adj) #Arrange by adjusted p.values
#Combine Kruskal-Wallis results at the family level for all days with sequence data:----
family_stats_dn1to9_combined <- rbind(kw_sig_dn1, kw_sig_d0, kw_sig_d1, kw_sig_d2, kw_sig_d3, kw_sig_d4,
                                      kw_sig_d5, kw_sig_d6, kw_sig_d7, kw_sig_d8, kw_sig_d9) %>% 
  filter(p.value.adj <= 0.05) %>% 
  write_tsv(path = paste0("data/process/family_kw_stats_dn1to9.tsv")) #Save combined dataframe as a .tsv
#Also write results to supplemental table excel file
write_xlsx(family_stats_dn1to9_combined, "submission/table_S17_family_kw_stats_dn1to9.xlsx", format_headers = FALSE)

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
    mutate(otu=factor(otu, levels=otus)) %>% 
    mutate(agg_rel_abund = agg_rel_abund + 1/10874) %>% # 10,874 is 2 times the subsampling parameter of 5437
    ggplot(aes(x= otu, y=agg_rel_abund, color=vendor))+
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_vendors,
                        labels=color_vendors)+
    geom_hline(yintercept=1/5437, color="gray")+
    geom_boxplot(outlier.shape = NA, size = 1.2)+
    #  geom_jitter(shape=19, size=2, alpha=0.6, position=position_jitterdodge(dodge.width=0.7, jitter.width=0.2)) + #Too busy when indiv. mice are shown
    labs(title=NULL, 
         x=NULL,
         y="Relative abundance (%)")+
    scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100), limits = c(1/10900, 1))+
    coord_flip()+
    theme_classic()+
    theme(plot.title=element_text(hjust=0.5),
          legend.position = "bottom",
          axis.text.y = element_text(face = "italic"), #Have the OTUs show up as italics
          text = element_text(size = 16)) # Change font size for entire plot
}

#Plots of the relative abundances of OTUs that significantly varied across sources of mice from day -1 to day 1----
Dn1toD1_otus_dn1 <- plot_otus_dx(shared_sig_otus_Dn1toD1, -1)+
  ggtitle("Baseline")+ #Title plot
  theme(plot.title = element_text(hjust = 0.5)) #Center plot title
save_plot("results/figures/Dn1toD1_otus_dn1.png", Dn1toD1_otus_dn1, base_height = 7, base_width = 8)
Dn1toD1_otus_d0 <- plot_otus_dx(shared_sig_otus_Dn1toD1, 0) +
  ggtitle("Clindamycin")+ #Title plot
  theme(plot.title = element_text(hjust = 0.5)) #Center plot title
save_plot("results/figures/Dn1toD1_otus_d0.png", Dn1toD1_otus_d0, base_height = 7, base_width = 8)
Dn1toD1_otus_d1 <- plot_otus_dx(shared_sig_otus_Dn1toD1, 1) +
  ggtitle("Post-infection")+ #Title plot
  theme(plot.title = element_text(hjust = 0.5)) #Center plot title
save_plot("results/figures/Dn1toD1_otus_d1.png", Dn1toD1_otus_d1, base_height = 7, base_width = 8)

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

#Combine OTU statistics for differences across souces of mice for days -1, 0 and 1:
otu_stats_dn1to1_combined <- rbind(kw_w_sig_dn1, kw_w_sig_d0, kw_w_sig_d1) %>% 
  write_tsv(path = paste0("data/process/otu_stats_dn1to1_combined.tsv")) %>%  #Save combined dataframe as a .tsv
  #Also write results to supplemental table excel file
  write_xlsx("submission/table_S8_otu_stats_dn1to1.xlsx", format_headers = FALSE)

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
#Also write results to supplemental table excel file
write_xlsx(otu_stats_dn1to9_combined, "submission/table_S16_otu_kw_stats_dn1to9.xlsx", format_headers = FALSE)

#Wilcoxan Signed rank test for relative abundance differences after clindamycin treatment (for all mice with paired data for day -1 versus day 0) at different taxonomic levels with Benjamini-Hochburg correction----
#Pull mice ids that have sequence data for both day -1 and day 0:
mice_seq_dn1_0_pairs <- agg_family_data %>% 
  filter(family == "Porphyromonadaceae") %>% #Random pick just to figure out what mice have sequence data
  filter(day == -1 | day == 0) %>% 
  filter(duplicated(mouse_id)) %>% #Pull mouse ids with sequence data for both day -1 and day 0
  pull(mouse_id) 

#Dataframe for statistical test at the family level
paired_family <- agg_family_data %>% 
  filter(mouse_id %in% mice_seq_dn1_0_pairs) %>% #Only select pairs with data for day -1 & day 0
  filter(day == -1 | day == 0) %>% #Experiment days that represent initial community and community post clindamycin treatment
  mutate(day = as.factor(day)) %>% 
  select(day, family, agg_rel_abund)

#Wilcoxon signed rank test for all day -1, day 0 pairs at the family level:
f_dn1to0_pairs <- paired_family %>% 
  group_by(family) %>% 
  nest() %>% 
  mutate(model=map(data, ~wilcox.test(.x$agg_rel_abund ~ .x$day, paired = TRUE) %>% tidy())) %>% 
  mutate(median = map(data, get_rel_abund_median_day)) %>% 
  unnest(c(model, median)) %>% 
  ungroup() 
#Adjust p-values for testing multiple families
f_dn1to0_pairs_stats_adjust <- f_dn1to0_pairs %>% 
  select(family, statistic, p.value, method, alternative, `-1`, `0`) %>% 
  mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
  arrange(p.value.adj) %>% 
  write_tsv(path = "data/process/family_stats_dn1to0.tsv")  
#Also write results to supplemental table excel file
write_xlsx(f_dn1to0_pairs_stats_adjust, "submission/table_S11_family_stats_dn1to0.xlsx", format_headers = FALSE)

#Make a list of significant families impacted by clindamycin treatment  
sig_family_pairs <- pull_significant_taxa(f_dn1to0_pairs_stats_adjust, family)
#18 families
sig_family_pairs_top10 <- sig_family_pairs[1:10]

#Plot of the families with significantly different relative abundances post clindamycin treatment across all sources of mice:----
clind_impacted_families_plot_dx <- function(timepoint){
  agg_family_data %>%
    filter(mouse_id %in% mice_seq_dn1_0_pairs) %>% #Only select pairs with data for day -1 & day 0
    filter(family %in% sig_family_pairs_top10) %>% 
    filter(day == timepoint) %>% 
    mutate(family=factor(family, levels=sig_family_pairs_top10)) %>% 
    mutate(agg_rel_abund = agg_rel_abund + 1/10874) %>% # 10,874 is 2 times the subsampling parameter of 5437
    ggplot(aes(x= family, y=agg_rel_abund, color=vendor))+
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_vendors,
                        labels=color_vendors)+
    geom_hline(yintercept=1/5437, color="gray")+
    geom_boxplot(outlier.shape = NA, size = 1.2)+
    #  geom_jitter(shape=19, size=2, alpha=0.6, position=position_jitterdodge(dodge.width=0.7, jitter.width=0.2)) +
    labs(title=NULL, 
         x=NULL,
         y="Relative abundance (%)")+
    scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100), limits = c(1/10900, 1))+
    coord_flip()+
    theme_classic()+
    theme(plot.title=element_text(hjust=0.5),
          text = element_text(size = 16),# Change font size for entire plot
          axis.text.y = element_text(face = "italic", size = 18), #Have the families show up as italics
          strip.background = element_blank(),
          legend.position = "bottom") 
}

clind_impacted_families_plot_dn1 <- clind_impacted_families_plot_dx(-1)+
  ggtitle("Baseline")+ #Title plot
  theme(plot.title = element_text(hjust = 0.5)) #Center plot title
save_plot(filename = paste0("results/figures/clind_impacted_families_plot_dn1.png"), clind_impacted_families_plot_dn1, base_height = 12, base_width = 7)
clind_impacted_families_plot_d0 <- clind_impacted_families_plot_dx(0)+
  ggtitle("Clindamycin")+ #Title plot
  theme(plot.title = element_text(hjust = 0.5)) #Center plot title
save_plot(filename = paste0("results/figures/clind_impacted_families_plot_d0.png"), clind_impacted_families_plot_d0, base_height = 12, base_width = 7)

#Dataframe for statistical test at the OTU level
paired_otu <- agg_otu_data %>% 
  filter(mouse_id %in% mice_seq_dn1_0_pairs) %>% #Only select pairs with data for day -1 & day 0
  filter(day == -1 | day == 0) %>% #Experiment days that represent initial community and community post clindamycin treatment
  mutate(day = as.factor(day)) %>% 
  select(day, otu, agg_rel_abund)

#Wilcoxon signed rank test for all day -1, day 0 pairs at the family level:
o_dn1to0_pairs <- paired_otu %>% 
  group_by(otu) %>% 
    nest() %>% 
    mutate(model=map(data, ~wilcox.test(.x$agg_rel_abund ~ .x$day, paired = TRUE) %>% tidy())) %>% 
    mutate(median = map(data, get_rel_abund_median_day)) %>% 
    unnest(c(model, median)) %>% 
    ungroup() 
#Adjust p-values for testing multiple families
o_dn1to0_pairs_stats_adjust <- o_dn1to0_pairs %>% 
    select(otu, statistic, p.value, method, alternative, `-1`, `0`) %>% 
    mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
    arrange(p.value.adj) %>% 
  write_tsv(path = "data/process/otu_stats_dn1to0.tsv") 
#Also write results to supplemental table excel file
write_xlsx(o_dn1to0_pairs_stats_adjust, "submission/table_S10_otu_stats_dn1to0.xlsx", format_headers = FALSE)

#Make a list of significant OTUs impacted by clindamycin treatment----  
sig_otu_pairs <- pull_significant_taxa(o_dn1to0_pairs_stats_adjust, otu)
# 153 OTUs
sig_otu_pairs_top10 <- sig_otu_pairs[1:10]

#Plot of the otus with significantly different relative abundances post clindamycin treatment across all sources of mice:----
clind_impacted_otus_plot_dx <- function(timepoint){ 
  agg_otu_data %>% 
  filter(mouse_id %in% mice_seq_dn1_0_pairs) %>% #Only select pairs with data for day -1 & day 0
  filter(otu %in% sig_otu_pairs_top10) %>% 
  filter(day == timepoint) %>% 
  mutate(otu=factor(otu, levels=sig_otu_pairs_top10)) %>% 
  mutate(agg_rel_abund = agg_rel_abund + 1/10874) %>% # 10,874 is 2 times the subsampling parameter of 5437
  ggplot(aes(x= otu, y=agg_rel_abund, color=vendor))+
  scale_colour_manual(name=NULL,
                      values=color_scheme,
                      breaks=color_vendors,
                      labels=color_vendors)+
  geom_hline(yintercept=1/5437, color="gray")+
  geom_boxplot(outlier.shape = NA, size = 1.2)+
  #  geom_jitter(shape=19, size=2, alpha=0.6, position=position_jitterdodge(dodge.width=0.7, jitter.width=0.2)) +
  labs(title=NULL, 
       x=NULL,
       y="Relative abundance (%)")+
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100), limits = c(1/10900, 1))+
  coord_flip()+
  theme_classic()+
  theme(plot.title=element_text(hjust=0.5),
        text = element_text(size = 16),# Change font size for entire plot
        axis.text.y = element_text(face = "italic", size = 18), #Have the families show up as italics
        strip.background = element_blank(),
        legend.position = "bottom") 
} 

clind_impacted_otus_plot_dn1 <- clind_impacted_otus_plot_dx(-1)+
  ggtitle("Baseline")+ #Title plot
  theme(plot.title = element_text(hjust = 0.5)) #Center plot title
save_plot(filename = paste0("results/figures/clind_impacted_otus_plot_dn1.png"), clind_impacted_otus_plot_dn1, base_height = 12, base_width = 7)
clind_impacted_otus_plot_d0 <- clind_impacted_otus_plot_dx(0)+
  ggtitle("Clindamycin")+ #Title plot
  theme(plot.title = element_text(hjust = 0.5)) #Center plot title
save_plot(filename = paste0("results/figures/clind_impacted_otus_plot_d0.png"), clind_impacted_otus_plot_d0, base_height = 12, base_width = 7)

#Comparison of Figure 4 (varied across colony source) and 5 (altered by clindamycin treatment) taxa----
Fig4_v_5_otus <- intersect_all(`shared_sig_otus_Dn1toD1`, `sig_otu_pairs_top10`)
#O OTUs overlap
Fig4_v_all_clind <- intersect_all(`shared_sig_otus_Dn1toD1`, `sig_otu_pairs`)
# 3 OTUs Lachnospiraceae (OTU 130), Lactobacillus (OTU 6), Enterococcus (OTU 23)
Fig4_v_5_families <- intersect_all(`shared_sig_families_Dn1toD1`, `sig_family_pairs_top10`)
#3 families overlap: "Porphyromonadaceae", "Enterococcaceae", "Lachnospiraceae"  
Fig4_v_all_clind <- intersect_all(`shared_sig_families_Dn1toD1`, `sig_family_pairs`)
#3 families overlap: "Porphyromonadaceae", "Enterococcaceae", "Lachnospiraceae"  

#Function to plot OTUs of interest that overlap with top 20 OTUs in 3 logistic regression models----
#Customize the x_annotation, y_position, and label argument values for each OTU prior to running the function
otu_over_time <- function(otu_plot, x_annotation, y_position, label){
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
    geom_point(otu_mice, mapping = aes(x=day, y=agg_rel_abund, color=vendor, shape=clearance_status_d7), alpha = .8, size  = 1.5, position = position_dodge(width = 0.6))+
    geom_line(otu_median, mapping = aes(x=day, y=median, color=vendor), size = 1, show.legend = FALSE)+
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_vendors,
                        labels=color_vendors)+
    scale_shape_manual(name="Cleared by Day 7",
                        values=c(4, 19),
                        breaks=c("colonized", "not_detectable"),
                        labels=c("no", "yes"), 
                        drop=FALSE, na.translate = TRUE, na.value = 1)+
    geom_hline(yintercept=1/5437, color="gray")+
    labs(title=otu_plot,
         x="Day",
         y="Relative abundance (%)") +
    scale_x_continuous(breaks = c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
                       limits = c(-1.5, 9.5)) +
    scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
    annotate("text", y = y_position, x = x_annotation, label = label, size =8)+
    theme_classic()+
    theme(plot.title=element_text(hjust=0.5, face = "italic"),
          text = element_text(size = 16)) # Change font size for entire plot
  save_plot(filename = paste0("results/figures/", otu_plot,"_time.png"), otu_time, base_aspect_ratio = 2)
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
otu_over_time(otu_plot = "Bacteroides (OTU 2)", x_annotation = x_OTU2, y_position = y_OTU2, label = label_OTU2)
#Set up statistical annotation arguments for Enterococcus (OTU 23):
x_OTU23 <- plot_otu_stats_dn1to9 %>% 
  filter(otu == "Enterococcus (OTU 23)") %>% pull(day)
y_OTU23 <- .4
label_OTU23 <- plot_otu_stats_dn1to9 %>% 
  filter(otu == "Enterococcus (OTU 23)") %>% pull(p.signif)
otu_over_time("Enterococcus (OTU 23)", x_annotation = x_OTU23, y_position = y_OTU23, label = label_OTU23)

#OTUs impacted by clindamycin treatment and are important features in at least 2/3 logistic regresssion models---- 
#Set up statistical annotation arguments for Enterobacteriaceae (OTU 1):
x_OTU1 <- plot_otu_stats_dn1to9 %>% 
  filter(otu == "Enterobacteriaceae (OTU 1)") %>% pull(day)
y_OTU1 <- 1.5
label_OTU1 <- plot_otu_stats_dn1to9 %>% 
  filter(otu == "Enterobacteriaceae (OTU 1)") %>% pull(p.signif)
otu_over_time("Enterobacteriaceae (OTU 1)", x_annotation = x_OTU1, y_position = y_OTU1, label = label_OTU1)
otu_over_time("Enterococcus (OTU 23)", x_annotation = x_OTU23, y_position = y_OTU23, label = label_OTU23) #Also varies across mouse sources
#Set up statistical annotation arguments for Porphyromonadaceae (OTU 7):
x_OTU7 <- plot_otu_stats_dn1to9 %>% 
  filter(otu == "Porphyromonadaceae (OTU 7)") %>% pull(day)
y_OTU7 <- .4
label_OTU7 <- plot_otu_stats_dn1to9 %>% 
  filter(otu == "Porphyromonadaceae (OTU 7)") %>% pull(p.signif)
otu_over_time("Porphyromonadaceae (OTU 7)", x_annotation = x_OTU7, y_position = y_OTU7, label = label_OTU7)

#Function to plot families of interest that overlap with top 20 OTUs in 3 logistic regression models 
#Customize the x_annoation, y_position, and label argument values for each family prior to running the function
family_over_time <- function(family_plot, x_annotation, y_position, label){
  family_median <- agg_family_data %>% 
    filter(family == family_plot) %>% 
    group_by(vendor, day) %>% 
    summarize(median=(median(agg_rel_abund + 1/10874))) %>% 
    ungroup
  family_mice <-  agg_family_data %>% 
    filter(family == family_plot) %>% 
    mutate(agg_rel_abund = agg_rel_abund + 1/10874) %>%
    select(vendor, day, agg_rel_abund, family, experiment, clearance_status_d7)
  family_time <- ggplot(NULL)+
    geom_point(family_mice, mapping = aes(x=day, y=agg_rel_abund, color=vendor, shape=clearance_status_d7), alpha = .8, size  = 1.5, position = position_dodge(width = 0.6))+
    geom_line(family_median, mapping = aes(x=day, y=median, color=vendor), size = 1, show.legend = FALSE)+
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_vendors,
                        labels=color_vendors)+
    scale_shape_manual(name="Cleared by Day 7",
                       values=c(4, 19),
                       breaks=c("colonized", "not_detectable"),
                       labels=c("no", "yes"), 
                       drop=FALSE, na.translate = TRUE, na.value = 1)+
    geom_hline(yintercept=1/5437, color="gray")+
    labs(title=family_plot,
         x="Day",
         y="Relative abundance (%)") +
    scale_x_continuous(breaks = c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
                       limits = c(-1.5, 9.5)) +
    scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
    annotate("text", y = y_position, x = x_annotation, label = label, size =8)+
    theme_classic()+
    theme(plot.title=element_text(hjust=0.5, face = "italic"),
          text = element_text(size = 16)) # Change font size for entire plot
  save_plot(filename = paste0("results/figures/", family_plot,"_time.png"), family_time, base_aspect_ratio = 2)
}

#Create column with statistical significance symbols for all OTUs:
plot_family_stats_dn1to9 <- family_stats_dn1to9_combined %>% 
  mutate(p.value.adj=round(p.value.adj, digits = 4)) %>% 
  filter(p.value.adj <= 0.05) %>%
  mutate(p.signif = case_when(
    p.value.adj > 0.05 ~ "NS",
    p.value.adj <= 0.05 ~ "*"
  )) 
#Families that differ across mouse sources at multiple timepoints and are important features in at least 2/3 logistic regresssion models---- 
#Set up statistical annotation arguments for "Bacteroidaceae":
x_Bacteroidaceae <- plot_family_stats_dn1to9 %>% 
  filter(family == "Bacteroidaceae") %>% pull(day)
y_Bacteroidaceae <- 1.5
label_Bacteroidaceae <- plot_family_stats_dn1to9 %>% 
  filter(family == "Bacteroidaceae") %>% pull(p.signif)
family_over_time("Bacteroidaceae", x_annotation = x_Bacteroidaceae, y_position = y_Bacteroidaceae, label = label_Bacteroidaceae)
#Set up statistical annotation arguments for "Bacteroidaceae":
x_Bacteroidaceae <- plot_family_stats_dn1to9 %>% 
  filter(family == "Bacteroidaceae") %>% pull(day)
y_Bacteroidaceae <- 1.5
label_Bacteroidaceae <- plot_family_stats_dn1to9 %>% 
  filter(family == "Bacteroidaceae") %>% pull(p.signif)
family_over_time("Bacteroidaceae", x_annotation = x_Bacteroidaceae, y_position = y_Bacteroidaceae, label = label_Bacteroidaceae)
#Set up statistical annotation arguments for "Deferribacteraceae":
x_Deferribacteraceae <- plot_family_stats_dn1to9 %>% 
  filter(family == "Deferribacteraceae") %>% pull(day)
y_Deferribacteraceae <- .4
label_Deferribacteraceae <- plot_family_stats_dn1to9 %>% 
  filter(family == "Deferribacteraceae") %>% pull(p.signif)
family_over_time("Deferribacteraceae", x_annotation = x_Deferribacteraceae, y_position = y_Deferribacteraceae, label = label_Deferribacteraceae)
#Set up statistical annotation arguments for "Enterococcaceae":
x_Enterococcaceae <- plot_family_stats_dn1to9 %>% 
  filter(family == "Enterococcaceae") %>% pull(day)
y_Enterococcaceae <- .4
label_Enterococcaceae <- plot_family_stats_dn1to9 %>% 
  filter(family == "Enterococcaceae") %>% pull(p.signif)
family_over_time("Enterococcaceae",  x_annotation = x_Enterococcaceae, y_position = y_Enterococcaceae, label = label_Enterococcaceae) 
#Set up statistical annotation arguments for "Lachnospiraceae":
x_Lachnospiraceae <- plot_family_stats_dn1to9 %>% 
  filter(family == "Lachnospiraceae") %>% pull(day)
y_Lachnospiraceae <- 1.5
label_Lachnospiraceae <- plot_family_stats_dn1to9 %>% 
  filter(family == "Lachnospiraceae") %>% pull(p.signif)
family_over_time("Lachnospiraceae",  x_annotation = x_Lachnospiraceae, y_position = y_Lachnospiraceae, label = label_Lachnospiraceae)

#Families impacted by clindamycin treatment and are important features in at least 2/3 logistic regresssion models---- 
#Set up statistical annotation arguments for "Bifidobacteriaceae":
x_Bifidobacteriaceae <- plot_family_stats_dn1to9 %>% 
  filter(family == "Bifidobacteriaceae") %>% pull(day)
y_Bifidobacteriaceae <- .4
label_Bifidobacteriaceae <- plot_family_stats_dn1to9 %>% 
  filter(family == "Bifidobacteriaceae") %>% pull(p.signif)
family_over_time("Bifidobacteriaceae", x_annotation = x_Bifidobacteriaceae, y_position = y_Bifidobacteriaceae, label = label_Bifidobacteriaceae)
#Set up statistical annotation arguments for "Coriobacteriaceae":
x_Coriobacteriaceae <- plot_family_stats_dn1to9 %>% 
  filter(family == "Coriobacteriaceae") %>% pull(day)
y_Coriobacteriaceae <- 0.01
label_Coriobacteriaceae <- plot_family_stats_dn1to9 %>% 
  filter(family == "Coriobacteriaceae") %>% pull(p.signif)
family_over_time("Coriobacteriaceae", x_annotation = x_Coriobacteriaceae, y_position = y_Coriobacteriaceae, label = label_Coriobacteriaceae)
family_over_time("Enterococcaceae",  x_annotation = x_Enterococcaceae, y_position = y_Enterococcaceae, label = label_Enterococcaceae) #Also varies across mouse sources
family_over_time("Lachnospiraceae",  x_annotation = x_Lachnospiraceae, y_position = y_Lachnospiraceae, label = label_Lachnospiraceae) #Also varies across mouse sources
#Set up statistical annotation arguments for "Ruminococcaceae":
x_Ruminococcaceae <- plot_family_stats_dn1to9 %>% 
  filter(family == "Ruminococcaceae") %>% pull(day)
y_Ruminococcaceae <- .4
label_Ruminococcaceae <- plot_family_stats_dn1to9 %>% 
  filter(family == "Ruminococcaceae") %>% pull(p.signif)
family_over_time("Ruminococcaceae", x_annotation = x_Ruminococcaceae, y_position = y_Ruminococcaceae, label = label_Ruminococcaceae)
#Set up statistical annotation arguments for "Verrucomicrobiaceae":
x_Verrucomicrobiaceae <- plot_family_stats_dn1to9 %>% 
  filter(family == "Verrucomicrobiaceae") %>% pull(day)
y_Verrucomicrobiaceae <- .6
label_Verrucomicrobiaceae <- plot_family_stats_dn1to9 %>% 
  filter(family == "Verrucomicrobiaceae") %>% pull(p.signif)
family_over_time("Verrucomicrobiaceae", x_annotation = x_Verrucomicrobiaceae, y_position = y_Verrucomicrobiaceae, label = label_Verrucomicrobiaceae)

#Do shared taxa associated with d7 cleared/colonized status?
#Function to test for differences in relative abundances at the family level according to day 7 colonization status:
w_d7status_f <- function(timepoint){
  family_stats <- agg_family_data %>% 
    filter(day == timepoint) %>% 
    filter(!is.na(clearance_status_d7)) %>% #Remove mice where clearance status day 7 was NA
    mutate(clearance_status_d7 = as.factor(clearance_status_d7)) %>% 
    select(clearance_status_d7, family, agg_rel_abund) %>% 
    group_by(family) %>% 
    nest() %>% 
    mutate(model=map(data, ~wilcox.test(.x$agg_rel_abund ~ .x$clearance_status_d7, paired = FALSE) %>% tidy())) %>% 
    mutate(mean = map(data, get_rel_abund_mean_d7status)) %>% 
    unnest(c(model, mean)) %>% 
    ungroup() 
  #Adjust p-values for testing multiple families
  family_stats_adjust <- family_stats %>% 
    select(family, statistic, p.value, method, alternative, colonized, not_detectable) %>% 
    mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
    arrange(p.value.adj) %>% 
    write_tsv(path = paste0("data/process/family_stats_d7status_d", timepoint, ".tsv"))
}

#Test only for days where input communities were used to create classification models (Day -1, 0, 1):
model_input_days <- c(-1, 0, 1)
for (d in model_input_days){
  w_d7status_f(d)
  #Make a list of significant families that differ according to day 7 C. diff status for a specific day  
  stats <- read_tsv(file = paste0("data/process/family_stats_d7status_d", d, ".tsv"))
  name <- paste("sig_family_status_", d, sep = "") 
  assign(name, pull_significant_taxa(stats, family))
}

`sig_family_status_-1` # 0 families
sig_family_status_0 # 0 families
sig_family_status_1 # 0 families

#Function to test for differences in relative abundances at the OTU level according to day 7 colonization status:
w_d7status_o <- function(timepoint){
  otu_stats <- agg_otu_data %>% 
    filter (day == timepoint) %>% 
    filter(!is.na(clearance_status_d7)) %>% #Remove mice where clearance status day 7 was NA
    mutate(clearance_status_d7 = as.factor(clearance_status_d7)) %>% 
    select(clearance_status_d7, otu, agg_rel_abund) %>% 
    group_by(otu) %>% 
    nest() %>% 
    mutate(model=map(data, ~wilcox.test(.x$agg_rel_abund ~ .x$clearance_status_d7, paired = FALSE) %>% tidy())) %>% 
    mutate(mean = map(data,  get_rel_abund_mean_d7status)) %>% 
    unnest(c(model, mean)) %>% 
    ungroup() 
  #Adjust p-values for testing multiple OTus
  otu_stats_adjust <- otu_stats %>% 
    select(otu, statistic, p.value, method, alternative, colonized, not_detectable) %>% 
    mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
    arrange(p.value.adj) %>% 
    write_tsv(path = paste0("data/process/otu_stats_d7status_d", timepoint, ".tsv"))
}

#Test only for days where input communities were used to create classification models (Day -1, 0, 1):
model_input_days <- c(-1, 0, 1)
for (d in model_input_days){
  w_d7status_o(d)
  #Make a list of significant OTUs that differ according to day 7 C. diff status for a specific day  
  stats <- read_tsv(file = paste0("data/process/otu_stats_d7status_d", d, ".tsv"))
  name <- paste("sig_otu_status_", d, sep = "") 
  assign(name, pull_significant_taxa(stats, otu))
}

`sig_otu_status_-1` # 1 OTU "Ruminococcaceae (OTU 467)"
sig_otu_status_0 # 0 OTUs
sig_otu_status_1 # 0 OTUs

#Set up statistical annotation arguments for Ruminococcaceae (OTU 467):
x_OTU467 <- plot_otu_stats_dn1to9 %>% 
  filter(otu == "Ruminococcaceae (OTU 467)") %>% pull(day) #No timepoints were significant
y_OTU467 <- .4
label_OTU467 <- plot_otu_stats_dn1to9 %>% 
  filter(otu == "Ruminococcaceae (OTU 467)") %>% pull(p.signif)
otu_over_time("Ruminococcaceae (OTU 467)", x_annotation = x_OTU467, y_position = y_OTU467, label = label_OTU467)

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
  #Adjust p-values for testing multiple families
  otu_stats_adjust <- otu_stats %>% 
    select(otu, statistic, p.value, method, alternative, `1`, `2`) %>% 
    mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
    arrange(p.value.adj) %>% 
    write_tsv(path = paste0("data/process/otu_stats_experiment_", vendor_name, ".tsv"))
}

#Test only for days where input communities were used to create classification models (Day -1, 0, 1):
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

#Examine top 20 taxa from logistic regression model that performed the best: Communities after clindamycin treatment at the OTU level----

#List of top 20 taxa from day 0 OTU logistic regression model
interp_otus_d0 <- interp_otus  %>% filter(model_input_day == 0) %>% pull(OTU)

#Split plotting of taxa up by relative abundance

#Function to plot specific OTUs
plot_interp_otus_d0 <- function(otu_name, otu_stats){
  d0_otu_model_plot <- agg_otu_data %>% 
    filter(day == 0) %>% 
    filter(otu == otu_name) %>% 
    mutate(agg_rel_abund = agg_rel_abund + 1/10874) %>% 
    group_by(vendor) %>% 
    mutate(median=(median(agg_rel_abund))) %>% #create a column of median values for each group
    ungroup() %>% 
    ggplot(aes(x=vendor, y =agg_rel_abund, colour= vendor))+
    scale_x_discrete(guide = guide_axis(n.dodge = 2))+
    geom_errorbar(aes(ymax = median, ymin = median), color = "gray50", size = 1)+ #Add lines to indicate the median for each group to the plot
    geom_jitter(aes(shape = clearance_status_d7), size=2, alpha=0.6, show.legend = FALSE) +
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_vendors,
                        labels=color_vendors)+
    scale_shape_manual(name="Cleared by Day 7",
                       values=c(4, 19),
                       breaks=c("colonized", "not_detectable"),
                       labels=c("no", "yes"), 
                       drop=FALSE, na.translate = TRUE, na.value = 1)+
    geom_hline(yintercept=1/5437, color="gray")+
    labs(title=otu_name,
         x=NULL,
         y="Relative abundance (%)") +
    scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
    theme_classic()+
    theme(plot.title=element_text(hjust=0.5, face = "italic"),
          text = element_text(size = 16)) +# Change font size for entire plot
    stat_pvalue_manual(data = otu_stats, label = "p.adj", y.position = "y.position", size = 6, bracket.size = .6) #Add manual p values to taxa that were also significantly different across vendors
    save_plot(filename = paste0("results/figures/d0_model_otu_", otu_name, ".png"), d0_otu_model_plot, base_aspect_ratio = 2)
}

#OTUs with high relative abundances in at least one group
#Check for significant pairwise comparisons:
otu1_stats <- otu_day0_stats %>% 
  filter(otu == "Enterobacteriaceae (OTU 1)") %>% 
  filter(p.adj <= 0.05) %>%  #Only show comparisons that were significant. p.value, which was adjusted < 0.05)
  mutate(p.adj="*") %>% #Just indicate whether statistically significant, exact p.adj values are in supplemental table
  mutate(y.position = c(.3, .9, 1.8, 2.1, 2.4, .6, 1.2, 1.5, .6))
plot_interp_otus_d0("Enterobacteriaceae (OTU 1)", otu1_stats)
otu2_stats <- otu_day0_stats %>% 
  filter(otu == "Bacteroides (OTU 2)") %>% 
  filter(p.adj <= 0.05) %>%  #Only show comparisons that were significant. p.value, which was adjusted < 0.05)
  mutate(p.adj="*") %>% #Just indicate whether statistically significant, exact p.adj values are in supplemental table
  mutate(y.position = (1:n())*.3)
plot_interp_otus_d0("Bacteroides (OTU 2)", otu2_stats)
otu16_stats <- otu_day0_stats %>% 
  filter(otu == "Proteus (OTU 16)") %>% 
  filter(p.adj <= 0.05) %>%  #Only show comparisons that were significant. p.value, which was adjusted < 0.05)
  mutate(p.adj="*") %>% #Just indicate whether statistically significant, exact p.adj values are in supplemental table
  mutate(y.position = (1:n())*.3)
plot_interp_otus_d0("Proteus (OTU 16)", otu16_stats)

# 1% or less relative abundances:
otu18_stats <- otu_day0_stats %>% 
  filter(otu == "Lactobacillus (OTU 18)") %>% 
  filter(p.adj <= 0.05) %>%  #Only show comparisons that were significant. p.value, which was adjusted < 0.05)
  mutate(p.adj="*") %>% #Just indicate whether statistically significant, exact p.adj values are in supplemental table
  mutate(y.position = 0) #No significant pairwise comparisons
plot_interp_otus_d0("Lactobacillus (OTU 18)", otu18_stats)
otu99_stats <- otu_day0_stats %>% 
  filter(otu == "Clostridium (OTU 99)") %>% 
  filter(p.adj <= 0.05) %>%  #Only show comparisons that were significant. p.value, which was adjusted < 0.05)
  mutate(p.adj="*") %>% #Just indicate whether statistically significant, exact p.adj values are in supplemental table
  mutate(y.position = 0) #No significant pairwise comparisons
plot_interp_otus_d0("Clostridium (OTU 99)", otu99_stats)
otu9_stats <- otu_day0_stats %>% 
  filter(otu == "Lachnospiraceae (OTU 9)") %>% 
  filter(p.adj <= 0.05) %>%  #Only show comparisons that were significant. p.value, which was adjusted < 0.05)
  mutate(p.adj="*") %>% #Just indicate whether statistically significant, exact p.adj values are in supplemental table
  mutate(y.position = 0) #No significant pairwise comparisons
plot_interp_otus_d0("Lachnospiraceae (OTU 9)", otu9_stats)

#Placeholder stats dataframe that is empty, will not plot stats for rest of taxa
otu_stats_placeholder <- otu9_stats
#Close to undetectable except a few mice:
plot_interp_otus_d0("Alishewanella (OTU 776)", otu_stats_placeholder)
plot_interp_otus_d0("Clostridium (OTU 226)", otu_stats_placeholder)
plot_interp_otus_d0("Eisenbergiella (OTU 164)", otu_stats_placeholder)
plot_interp_otus_d0("Erysipelotrichaceae (OTU 234)", otu_stats_placeholder)
plot_interp_otus_d0("Lachnospiraceae (OTU 33)", otu_stats_placeholder)
plot_interp_otus_d0("Lachnospiraceae (OTU 38)", otu_stats_placeholder)
plot_interp_otus_d0("Lachnospiraceae (OTU 56)", otu_stats_placeholder)
plot_interp_otus_d0("Lactobacillus (OTU 834)", otu_stats_placeholder)
plot_interp_otus_d0("Porphyromonadaceae (OTU 7)", otu_stats_placeholder)
plot_interp_otus_d0("Porphyromonadaceae (OTU 22)", otu_stats_placeholder)
plot_interp_otus_d0("Porphyromonadaceae (OTU 54)", otu_stats_placeholder)
plot_interp_otus_d0("Ruminococcaceae (OTU 60)", otu_stats_placeholder)
plot_interp_otus_d0("Ruminococcaceae (OTU 520)", otu_stats_placeholder)
plot_interp_otus_d0("Escherichia/Shigella (OTU 610)", otu_stats_placeholder)# Have to do this one separately because of the slash
d0_escherichia_model_plot <- agg_otu_data %>% 
  filter(day == 0) %>% 
  filter(otu == "Escherichia/Shigella (OTU 610)") %>% 
  mutate(agg_rel_abund = agg_rel_abund + 1/10874) %>% 
  group_by(vendor) %>% 
  mutate(median=(median(agg_rel_abund))) %>% #create a column of median values for each group
  ungroup() %>% 
  ggplot(aes(x=vendor, y =agg_rel_abund, colour= vendor))+
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  geom_errorbar(aes(ymax = median, ymin = median), color = "gray50", size = 1)+ #Add lines to indicate the median for each group to the plot
  geom_jitter(aes(shape = clearance_status_d7), size=2, alpha=0.6, show.legend = FALSE) +
  scale_colour_manual(name=NULL,
                      values=color_scheme,
                      breaks=color_vendors,
                      labels=color_vendors)+
  scale_shape_manual(name="Cleared by Day 7",
                     values=c(4, 19),
                     breaks=c("colonized", "not_detectable"),
                     labels=c("no", "yes"), 
                     drop=FALSE, na.translate = TRUE, na.value = 1)+
  geom_hline(yintercept=1/5437, color="gray")+
  labs(title="Escherichia/Shigella (OTU 610)",
       x=NULL,
       y="Relative abundance (%)") +
  scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
  theme_classic()+
  theme(plot.title=element_text(hjust=0.5, face = "italic"),
        text = element_text(size = 16)) # Change font size for entire plot
save_plot(filename = paste0("results/figures/d0_model_otu_escherichia-shigella.png"), d0_escherichia_model_plot, base_aspect_ratio = 2)

#Plot of all 20 OTus together using facet_wrap
d0_otu_model_taxa <- agg_otu_data %>% 
    filter(day == 0) %>% 
    filter(otu %in% interp_otus_d0) %>% 
    mutate(agg_rel_abund = agg_rel_abund + 1/10874) %>% 
    group_by(vendor, otu) %>% 
    mutate(median=(median(agg_rel_abund))) %>% #create a column of median values for each group
    ungroup() %>% 
    ggplot(aes(x=vendor, y =agg_rel_abund, colour= vendor))+
    scale_x_discrete(guide = guide_axis(n.dodge = 2))+
    geom_errorbar(aes(ymax = median, ymin = median), color = "gray50", size = 1)+ #Add lines to indicate the median for each group to the plot
    geom_jitter(aes(shape = experiment), size=2, alpha=0.6, show.legend = FALSE) +
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_vendors,
                        labels=color_vendors)+
    scale_shape_manual(name="Cleared by Day 7",
                       values=c(4, 19),
                       breaks=c("colonized", "not_detectable"),
                       labels=c("no", "yes"), 
                       drop=FALSE, na.translate = TRUE, na.value = 1)+
    geom_hline(yintercept=1/5437, color="gray")+
    facet_wrap(~ otu)


