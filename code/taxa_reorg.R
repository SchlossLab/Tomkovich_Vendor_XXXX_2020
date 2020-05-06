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
                                              "_unclassified" = " Unclassified"))) %>% 
  # Separate taxonomic levels into separate columns according to semi-colon.
  separate(taxonomy, into=c("kingdom", "phylum", "class", "order", "family", "genus"), sep=';')

# Import otu_data for samples
otu_data <- read_tsv("data/process/vendors.subsample.shared", col_types=cols(Group=col_character())) %>% 
  select(-label, -numOtus) %>% 
  rename(id=Group) %>% 
  gather(-id, key="otu", value="count") %>% 
  mutate(rel_abund=count/5437) #Use 5437, because this is the subsampling parameter chosen.)

#Merge otu_data to taxonomy data frame
agg_taxa_data <- inner_join(otu_data, taxonomy)

# Function to summarize relative abundance level for a given taxonomic level (ex. genus, family, phlyum, etc.)
agg_taxonomic_data <- function(taxonomic_level) {
  taxonomic_level <- enquo(taxonomic_level) #Part of transformation of taxonomic_level argument into a column name
  agg_taxa_data %>% 
    group_by(id, !!taxonomic_level) %>% #!!Completes the transformation of taxonomic_level argument into a column name
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

# Relative abundance data at the genus level:
agg_genus_data <- agg_taxonomic_data(genus)
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

#Function to find which significant otus/genera/families are shared 
intersect_all <- function(a,b,...){
  Reduce(intersect, list(a,b,...))
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
  mutate(mean = map(data, get_rel_abund_mean_vendor)) %>% 
  unnest(c(model, mean)) %>% 
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

#Shared significant familes across days, but excluding day 2 and day 8 (low number of samples sequenced)----
shared_sig_families <- intersect_all(`sig_family_day-1`, sig_family_day0, sig_family_day1, sig_family_day3, sig_family_day4, sig_family_day5, sig_family_day6, sig_family_day7, sig_family_day9)
# 6 families: Betaproteobacteria Unclassified, Burkholderiales Unclassified, Sutterellaceae, Deferribacteraceae, Bacteroidaceae, Porphyromonadaceae

#Shared significant genera across D-1 to D1----
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
        axis.text.y = element_text(face = "italic"), #Have the families show up as italics
        text = element_text(size = 16)) # Change font size for entire plot
}

#Plots of the relative abundances of families that significantly vary across sources of mice from day -1 to day 1----
Dn1toD1_families_dn1 <- plot_families_dx(shared_sig_families_Dn1toD1, -1) 
save_plot("results/figures/Dn1toD1_families_dn1.png", Dn1toD1_families_dn1, base_aspect_ratio = 2)
Dn1toD1_families_d0 <- plot_families_dx(shared_sig_families_Dn1toD1, 0) 
save_plot("results/figures/Dn1toD1_families_d0.png", Dn1toD1_families_d0, base_aspect_ratio = 2)
Dn1toD1_families_d1 <- plot_families_dx(shared_sig_families_Dn1toD1, 1) 
save_plot("results/figures/Dn1toD1_families_d1.png", Dn1toD1_families_d1, base_aspect_ratio = 2)

# Perform pairwise Wilcoxan rank sum tests for families that were significantly different across sources of mice on a series of days----
pairwise_day_family <- function(timepoint, sig_family_dayX){
  family_stats <- agg_family_data %>% 
    filter(day == timepoint) %>%
    select(vendor, family, agg_rel_abund) %>% 
    group_by(family) %>% 
    nest() %>% 
    mutate(model=map(data, ~kruskal.test(x=.x$agg_rel_abund, g=as.factor(.x$vendor)) %>% tidy())) %>% 
    mutate(mean = map(data, get_rel_abund_mean_vendor)) %>% 
    unnest(c(model, mean)) %>% 
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

#Combine OTU statistics for differences across souces of mice for days -1, 0 and 1:
family_stats_dn1to1_combined <- rbind(kw_w_sig_dn1, kw_w_sig_d0, kw_w_sig_d1) %>% 
  write_tsv(path = paste0("data/process/family_stats_dn1to1_combined.tsv")) #Save combined dataframe as a .tsv

#Function to test at the genus level:
kruskal_wallis_g <- function(timepoint){
  genus_stats <- agg_genus_data %>% 
    filter(day == timepoint) %>%
    select(vendor, genus, agg_rel_abund) %>% 
    group_by(genus) %>% 
    nest() %>% 
    mutate(model=map(data, ~kruskal.test(x=.x$agg_rel_abund, g=as.factor(.x$vendor)) %>% tidy())) %>% 
    mutate(mean = map(data, get_rel_abund_mean_vendor)) %>% 
    unnest(c(model, mean)) %>% 
    ungroup() 
  #Adjust p-values for testing multiple genera
  genus_stats_adjust <- genus_stats %>% 
    select(genus, statistic, p.value, parameter, method, Schloss, Young, Jackson, `Charles River`, Taconic, Envigo) %>% 
    mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
    arrange(p.value.adj) %>% 
    write_tsv(path = paste0("data/process/genus_stats_day_", timepoint, ".tsv"))
}
# Perform kruskal wallis tests at the genus level for all days of the experiment that were sequenced----
for (d in exp_days_seq){
  kruskal_wallis_g(d)
  #Make a list of significant genera across sources of mice for a specific day  
  stats <- read_tsv(file = paste0("data/process/genus_stats_day_", d, ".tsv"))
  name <- paste("sig_genus_day", d, sep = "") 
  assign(name, pull_significant_taxa(stats, genus))
}

#Shared significant genera across all days except 2 and 8 (low number of samples sequenced that day)----
shared_sig_genera <- intersect_all(`sig_genus_day-1`, sig_genus_day0, sig_genus_day1, sig_genus_day3, sig_genus_day4, sig_genus_day5, sig_genus_day6, sig_genus_day7, sig_genus_day9)
# 9 genera: Betaproteobacteria Unclassified, Burkholderiales Unclassified, Parasutterella, Parabacteroides, Mucispirillum, Turicibacter, Bacteroides, Proteus, Clostridium XVIII

#Shared significant genera across D-1 to D1----
shared_sig_genera_Dn1toD1 <- intersect_all(`sig_genus_day-1`, sig_genus_day0, sig_genus_day1)
# 11 genera: Betaproteobacteria Unclassified, Burkholderiales Unclassified, Parasutterella, Parabacteroides, Mucispirillum, Turicibacter, Bacteroides, Clostridium XVIII, Enterococcus, Lachnospiraceae Unclassified

#Function to test at the otu level:
kruskal_wallis_otu <- function(timepoint){
  otu_stats <- agg_otu_data %>% 
    filter(day == timepoint) %>%
    select(vendor, otu, agg_rel_abund) %>% 
    group_by(otu) %>% 
    nest() %>% 
    mutate(model=map(data, ~kruskal.test(x=.x$agg_rel_abund, g=as.factor(.x$vendor)) %>% tidy())) %>% 
    mutate(mean = map(data, get_rel_abund_mean_vendor)) %>% 
    unnest(c(model, mean)) %>% 
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

#Shared significant OTUs across all days except 2 and 8 (low number of samples sequenced that day)----
shared_sig_otus <- intersect_all(`sig_otu_day-1`, sig_otu_day0, sig_otu_day1, sig_otu_day3, sig_otu_day4, sig_otu_day5, sig_otu_day6, sig_otu_day7, sig_otu_day9)
# 9 OTUs: [1] "Parasutterella (OTU 26)" "Burkholderiales (OTU 34)" "Betaproteobacteria (OTU 58)" "Lactobacillus (OTU 49)"     
# "Parabacteroides (OTU 5)" "Lactobacillus (OTU 31)" "Bacteroides (OTU 2)" "Proteus (OTU 16)" "Ruminococcaceae (OTU 152)" 

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
          axis.text.y = element_text(face = "italic"), #Have the OTUs show up as italics
          text = element_text(size = 16)) # Change font size for entire plot
}

#Plots of the relative abundances of OTUs that significantly vary across sources of mice from day -1 to day 1----
Dn1toD1_otus_dn1 <- plot_otus_dx(shared_sig_otus_Dn1toD1, -1) 
save_plot("results/figures/Dn1toD1_otus_dn1.png", Dn1toD1_otus_dn1, base_height = 6, base_width = 8)
Dn1toD1_otus_d0 <- plot_otus_dx(shared_sig_otus_Dn1toD1, 0) 
save_plot("results/figures/Dn1toD1_otus_d0.png", Dn1toD1_otus_d0, base_height = 6, base_width = 8)
Dn1toD1_otus_d1 <- plot_otus_dx(shared_sig_otus_Dn1toD1, 1) 
save_plot("results/figures/Dn1toD1_otus_d1.png", Dn1toD1_otus_d1, base_height = 6, base_width = 8)

# Perform pairwise Wilcoxan rank sum tests for otus that were significantly different across sources of mice on a series of days----
pairwise_day_otu <- function(timepoint, sig_otu_dayX){
  otu_stats <- agg_otu_data %>% 
    filter(day == timepoint) %>%
    select(vendor, otu, agg_rel_abund) %>% 
    group_by(otu) %>% 
    nest() %>% 
    mutate(model=map(data, ~kruskal.test(x=.x$agg_rel_abund, g=as.factor(.x$vendor)) %>% tidy())) %>% 
    mutate(mean = map(data, get_rel_abund_mean_vendor)) %>% 
    unnest(c(model, mean)) %>% 
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
  write_tsv(path = paste0("data/process/otu_stats_dn1to1_combined.tsv")) #Save combined dataframe as a .tsv

# Plot OTUs of interest (overlap with top 20 otus that came out of logistic regression model built from corresponding input day community: -1, 0, or 1)----
#Function to plot 1 significant OTU relative abundance across sources of mice at a specific timepoint----
#Arguments:
# name = name of OTU to plot.
# timepoint = timepoint to be analyzed
plot_otu_timepoint <- function(name, timepoint, stats){
  plot_otu <- agg_otu_data %>% 
    filter(otu == name) %>% 
    filter(day == timepoint) %>% 
    mutate(otu=factor(otu, name)) %>% 
    mutate(agg_rel_abund = agg_rel_abund + 1/10874) %>% # 10,874 is 2 times the subsampling parameter of 5437
    ggplot(aes(x= vendor, y=agg_rel_abund, color=vendor))+
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_vendors,
                        labels=color_vendors)+
    geom_hline(yintercept=1/5437, color="gray")+
    geom_boxplot(outlier.shape = NA, size = 1.2)+
    geom_jitter(shape=19, size=2, alpha=0.6, position=position_jitterdodge(dodge.width=0.7, jitter.width=0.2)) +
    labs(title=name, 
         x=NULL,
         y="Relative abundance (%)")+
    scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
    stat_pvalue_manual(data = stats, label = "p.adj", y.position = "y.position")+
    theme_classic()+
    theme(plot.title=element_text(hjust=0.5, face = "italic"))+
    theme(legend.position = "none") + #Get rid of legend
    theme(text = element_text(size = 16))  # Change font size for entire plot
  save_plot(filename = paste0("exploratory/notebook/day", timepoint, "/", name, "_at_day", timepoint,".png"), plot_otu, base_height = 11, base_width = 8.5, base_aspect_ratio = 2)
}
#Format day -1 stats dataframe for Otu 189:
Otu189_dn1_stats <- otu_dayn1_stats %>% 
  filter(otu == "Erysipelotrichaceae (OTU 189)" & p.adj <= 0.05) %>% 
  mutate(p.adj=round(p.adj, digits = 4)) %>% #No significant pairwise comparisons
  mutate(y.position = NA)
#Otu 189 plot----
plot_otu_timepoint("Erysipelotrichaceae (OTU 189)", -1, Otu189_dn1_stats)

#Format day -1 stats dataframe for Otu 58:
Otu58_dn1_stats <- otu_dayn1_stats %>% 
  filter(otu == "Betaproteobacteria (OTU 58)" & p.adj <= 0.05) %>% 
  mutate(p.adj=round(p.adj, digits = 4)) %>% #No significant pairwise comparisons
  mutate(y.position = c(1, .8, .4, .2, 1e-20))
#Otu 58 plot----
plot_otu_timepoint("Betaproteobacteria (OTU 58)", -1, Otu58_dn1_stats)

#Format day -1 stats dataframe for Otu 23:
Otu23_dn1_stats <- otu_dayn1_stats %>% 
  filter(otu == "Enterococcus (OTU 23)" & p.adj <= 0.05) %>% 
  mutate(p.adj=round(p.adj, digits = 4)) %>% #No significant pairwise comparisons
  mutate(y.position = NA)
#Otu 23 plot----
plot_otu_timepoint("Enterococcus (OTU 23)", -1, Otu23_dn1_stats)

#Format day -1 stats dataframe for Otu 34:
Otu34_dn1_stats <- otu_dayn1_stats %>% 
  filter(otu == "Burkholderiales (OTU 34)" & p.adj <= 0.05) %>% 
  mutate(p.adj=round(p.adj, digits = 4)) %>% #No significant pairwise comparisons
  mutate(y.position = c(1, .8, .4, .2, 1e-20))
#Otu 34 plot----
plot_otu_timepoint("Burkholderiales (OTU 34)", -1, Otu34_dn1_stats)

#Format day -1 stats dataframe for Otu 293:
Otu293_dn1_stats <- otu_dayn1_stats %>% 
  filter(otu == "Coriobacteriaceae (OTU 293)" & p.adj <= 0.05) %>% 
  mutate(p.adj=round(p.adj, digits = 4)) %>% #No significant pairwise comparisons
  mutate(y.position = c(1.2, 1, .8, .4, .2, 1e-20))
#Otu 293 plot----
plot_otu_timepoint("Coriobacteriaceae (OTU 293)", -1, Otu293_dn1_stats)

#Format day 0 stats dataframe for Otu 1:
Otu1_d0_stats <- otu_day0_stats %>% 
  filter(otu == "Enterobacteriaceae (OTU 1)" & p.adj <= 0.05) %>% 
  mutate(p.adj=round(p.adj, digits = 4)) %>% 
  mutate(y.position = c(.6, .8, 1.7, 2, 2.2, .4, 1.1, 1.4, .1)) #9 pairwise sig.
#Otu 1 plot----
plot_otu_timepoint("Enterobacteriaceae (OTU 1)", 0, Otu1_d0_stats)

#Format day 0 stats dataframe for Otu 2:
Otu2_d0_stats <- otu_day0_stats %>% 
  filter(otu == "Bacteroides (OTU 2)" & p.adj <= 0.05) %>% 
  mutate(p.adj=round(p.adj, digits = 4)) %>% 
  mutate(y.position = c(.2, .6, .9, 1.2, 1.5)) #5 pairwise sig.
#Otu 2 plot----
plot_otu_timepoint("Bacteroides (OTU 2)", 0, Otu2_d0_stats)

#Format day 0 stats dataframe for Otu 16:
Otu16_d0_stats <- otu_day0_stats %>% 
  filter(otu == "Proteus (OTU 16)" & p.adj <= 0.05) %>% 
  mutate(p.adj=round(p.adj, digits = 4)) %>% 
  mutate(y.position = c(1.4, 1.1, .8, .5, .2)) #5 pairwise sig.
#Otu 16 plot----
plot_otu_timepoint("Proteus (OTU 16)", 0, Otu16_d0_stats)


#Wilcoxan Rank Sum test for relative abundance differences within a single source of mice after clindamycin treatment (day -1 versus day 0) at different taxonomic levels with Benjamini-Hochburg correction----
mouse_sources <- levels(metadata$vendor)
#Function to test at the family level:
w_day_f <- function(source){
  family_stats <- agg_family_data %>% 
    filter(vendor == source) %>%
    filter (day == -1 | day == 0) %>% #Experiment days that represent initial community and community post clindamycin treatment
    select(day, family, agg_rel_abund) %>% 
    group_by(family) %>% 
    nest() %>% 
    mutate(model=map(data, ~wilcox.test(x=.x$agg_rel_abund, g=as.factor(.x$day)) %>% tidy())) %>% 
    mutate(mean = map(data, get_rel_abund_mean_day)) %>% 
    unnest(c(model, mean)) %>% 
    ungroup() 
  #Adjust p-values for testing multiple families
  family_stats_adjust <- family_stats %>% 
    select(family, statistic, p.value, method, alternative, `-1`, `0`) %>% 
    mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
    arrange(p.value.adj) %>% 
    write_tsv(path = paste0("data/process/family_stats_dn1to0_", source, ".tsv"))
}

# Perform wilcoxan rank sum tests to test impact of clindamycin on relative abundances at the family level within each source of mice----
for (s in mouse_sources){
  w_day_f(s)
  #Make a list of significant genera across time for a specific source of mice  
  stats <- read_tsv(file = paste0("data/process/family_stats_dn1to0_", s, ".tsv"))
  name <- paste("sig_family_", s, sep = "") 
  assign(name, pull_significant_taxa(stats, family))
}

#Combine stats for each individual source for families that varied from day -1 to 1 with adjusted p.values < 0.05----
#Read in individual tables:
for (s in mouse_sources){
  name <- paste("w_sig_family_", s, sep = "") 
  assign(name, read_tsv(file = paste0("data/process/family_stats_dn1to0_", s, ".tsv")))
}
#Filter tables so that only significant taxa are shown and make a new column to note source
w_sig_family_Schloss <- w_sig_family_Schloss %>% 
  filter(p.value.adj < 0.05) %>% 
  mutate(source = "Schloss")
w_sig_family_Young <- w_sig_family_Young %>% 
  filter(p.value.adj < 0.05) %>% 
  mutate(source = "Young")
w_sig_family_Jackson <- w_sig_family_Jackson %>% 
  filter(p.value.adj < 0.05) %>% 
  mutate(source = "Jackson")
`w_sig_family_Charles River` <- `w_sig_family_Charles River` %>% 
  filter(p.value.adj < 0.05) %>% 
  mutate(source = "Charles River")
w_sig_family_Taconic <- w_sig_family_Taconic %>% 
  filter(p.value.adj < 0.05) %>% 
  mutate(source = "Taconic")
w_sig_family_Envigo <- w_sig_family_Envigo %>% 
  filter(p.value.adj < 0.05) %>% 
  mutate(source = "Envigo")
#Combine tables for all sources of mice
family_stats_sources_combined <- rbind(w_sig_family_Schloss, w_sig_family_Young, w_sig_family_Jackson, 
                                    `w_sig_family_Charles River`, w_sig_family_Taconic, w_sig_family_Envigo) %>% 
  arrange(family) %>% 
  write_tsv(path = paste0("data/process/family_dn1to0_stats_all_sources.tsv")) #Save combined dataframe as a .tsv

#Families with relative abundances significantly altered by clindamycin treatment that are shared across sources of mice----
#Shared families across all sources of mice:
shared_all_sources_families <- intersect_all(`sig_family_Schloss`, `sig_family_Young`, `sig_family_Jackson`, `sig_family_Charles River`, `sig_family_Taconic`, `sig_family_Envigo`)
summary(shared_all_sources_families) # 6 familes 
#"Lactobacillaceae", "Bacteroidaceae", "Enterobacteriaceae", "Lachnospiraceae", "Ruminococcaceae", "Porphyromonadaceae"

#Compare to the list of families with relative abundances that significantly vary across sources of mice on day -1, 0, and 1 from above:
families_shared <- intersect_all(`shared_all_sources_families`, `shared_sig_families_Dn1toD1`)
# "Bacteroidaceae", "Lachnospiraceae", "Porphyromonadaceae"

#Shared families across Schloss and Young lab mice:
shared_Schloss_Young_families <- intersect_all(`sig_family_Schloss`, `sig_family_Young`)
summary(shared_Schloss_Young_families) #14 families 

#Shared families across Schloss, Young and Charles River mice:
shared_Schloss_Young_CR_families <- intersect_all(`sig_family_Schloss`, `sig_family_Young`, `sig_family_Charles River`)
summary(shared_Schloss_Young_CR_families) #8 families 

# Shared families across Jackson, Taconic and Envigo mice:
shared_JAX_Tac_Env_families <- intersect_all(`sig_family_Jackson`, `sig_family_Taconic`, `sig_family_Envigo`)
summary(shared_JAX_Tac_Env_families) #6 families 

#Shared families across mice purchased from 4 vendors:
shared_4_vendors_families <- intersect_all(`sig_family_Jackson`, `sig_family_Charles River`, `sig_family_Taconic`, `sig_family_Envigo`)
summary(shared_4_vendors_families) #6 families 

#Shared families across Jackson and Charles River mice:
shared_JAX_CR_families <- intersect_all(`sig_family_Jackson`, `sig_family_Charles River`)
summary(shared_JAX_CR_families) #8 families

#Shared families across Taconic and Envigo mice:
shared_Tac_Env_families <- intersect_all(`sig_family_Taconic`, `sig_family_Envigo`)
summary(shared_Tac_Env_families) #7 families

#Plot of the families with significantly different relative abundances post clindamycin treatment across all sources of mice:----
clind_impacted_families_plot <- agg_family_data %>% 
    filter(family %in% shared_all_sources_families) %>% 
    filter(day == -1 | day == 0) %>% 
    mutate(family=factor(family, levels=shared_all_sources_families)) %>% 
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
    facet_wrap(~day, dir = "v")+
    theme(plot.title=element_text(hjust=0.5),
          axis.text.y = element_text(face = "italic"), #Have the families show up as italics
          text = element_text(size = 16),# Change font size for entire plot
          strip.background = element_blank(),
          legend.position = "bottom") 
save_plot(filename = paste0("results/figures/clind_impacted_families_plot.png"), clind_impacted_families_plot, base_height = 8, base_width = 5)
  

#Function to test at the genus level:
w_day_g <- function(source){
  genus_stats <- agg_genus_data %>% 
    filter(vendor == source) %>%
    filter (day == -1 | day == 0) %>% #Experiment days that represent initial community and community post clindamycin treatment
    select(day, genus, agg_rel_abund) %>% 
    group_by(genus) %>% 
    nest() %>% 
    mutate(model=map(data, ~wilcox.test(x=.x$agg_rel_abund, g=as.factor(.x$day)) %>% tidy())) %>% 
    mutate(mean = map(data, get_rel_abund_mean_day)) %>% 
    unnest(c(model, mean)) %>% 
    ungroup() 
  #Adjust p-values for testing multiple genera
  genus_stats_adjust <- genus_stats %>% 
    select(genus, statistic, p.value, method, alternative, `-1`, `0`) %>% 
    mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
    arrange(p.value.adj) %>% 
    write_tsv(path = paste0("data/process/genus_stats_dn1to0_", source, ".tsv"))
}

# Perform wilcoxan rank sum tests to test impact of clindamycin on relative abundances at the genus level within each source of mice----
for (s in mouse_sources){
  w_day_g(s)
  #Make a list of significant genera across time for a specific source of mice  
  stats <- read_tsv(file = paste0("data/process/genus_stats_dn1to0_", s, ".tsv"))
  name <- paste("sig_genus_", s, sep = "") 
  assign(name, pull_significant_taxa(stats, genus))
}

# Number of significant genera for each source of mice:
summary(`sig_genus_Schloss`) # 9 significant genera
summary(`sig_genus_Young`) # 22 significant genera
summary(`sig_genus_Jackson`) # 7 significant genera
summary(`sig_genus_Charles River`) # 12 significant genera
summary(`sig_genus_Taconic`) # 25 significant genera
summary(`sig_genus_Envigo`) # 12 significant genera

#Genera with relative abundances significantly altered by clindamycin treatment that are shared across sources of mice----
#Shared genera across all sources of mice:
shared_all_sources <- intersect_all(`sig_genus_Schloss`, `sig_genus_Young`, `sig_genus_Jackson`, `sig_genus_Charles River`, `sig_genus_Taconic`, `sig_genus_Envigo`)
summary(shared_all_sources) # 3 genera that significantly change over time and shared by all sources of mice
#"Lactobacillus", "Bacteroides", "Lachnospiraceae unclassified"

#Shared genera across Schloss and Young lab mice:
shared_Schloss_Young <- intersect_all(`sig_genus_Schloss`, `sig_genus_Young`)
summary(shared_Schloss_Young) #9 genera 

#Shared genera across Schloss, Young and Charles River mice:
shared_Schloss_Young_CR <- intersect_all(`sig_genus_Schloss`, `sig_genus_Young`, `sig_genus_Charles River`)
summary(shared_Schloss_Young_CR) #6 genera 

#Shared genera across Jackson, Taconic and Envigo mice:
shared_JAX_Tac_Env <- intersect_all(`sig_genus_Jackson`, `sig_genus_Taconic`, `sig_genus_Envigo`)
summary(shared_JAX_Tac_Env) #4 genera

#Shared genera across mice purchased from 4 vendors:
shared_4_vendors <- intersect_all(`sig_genus_Jackson`, `sig_genus_Charles River`, `sig_genus_Taconic`, `sig_genus_Envigo`)
summary(shared_4_vendors) #3 genera

#Shared genera across Jackson and Charles River mice:
shared_JAX_CR <- intersect_all(`sig_genus_Jackson`, `sig_genus_Charles River`)
summary(shared_JAX_CR) #4 genera

#Shared genera across Taconic and Envigo mice:
shared_Tac_Env <- intersect_all(`sig_genus_Taconic`, `sig_genus_Envigo`)
summary(shared_Tac_Env) #9 genera

#Function to test at the OTU level:
w_day_o <- function(source){
  otu_stats <- agg_otu_data %>% 
    filter(vendor == source) %>%
    filter (day == -1 | day == 0) %>% #Experiment days that represent initial community and community post clindamycin treatment
    select(day, otu, agg_rel_abund) %>% 
    group_by(otu) %>% 
    nest() %>% 
    mutate(model=map(data, ~wilcox.test(x=.x$agg_rel_abund, g=as.factor(.x$day)) %>% tidy())) %>% 
    mutate(mean = map(data, get_rel_abund_mean_day)) %>% 
    unnest(c(model, mean)) %>% 
    ungroup() 
  #Adjust p-values for testing multiple OTUs
  otu_stats_adjust <- otu_stats %>% 
    select(otu, statistic, p.value, method, alternative, `-1`, `0`) %>% 
    mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
    arrange(p.value.adj) %>% 
    write_tsv(path = paste0("data/process/otu_stats_dn1to0_", source, ".tsv"))
}

# Perform wilcoxan rank sum tests to test impact of clindamycin on relative abundances at the OTU level within each source of mice----
for (s in mouse_sources){
  w_day_o(s)
  #Make a list of significant OTUs across time for a specific source of mice  
  stats <- read_tsv(file = paste0("data/process/otu_stats_dn1to0_", s, ".tsv"))
  name <- paste("sig_otu_", s, sep = "") 
  assign(name, pull_significant_taxa(stats, otu))
}

# Number of significant genera for each source of mice:
summary(`sig_otu_Schloss`) # 2 significant OTUs
summary(`sig_otu_Young`) # 2 significant OTUs
summary(`sig_otu_Jackson`) # 0 significant OTUs
summary(`sig_otu_Charles River`) # 0 significant OTUs
summary(`sig_otu_Taconic`) # 0 significant OTUs
summary(`sig_otu_Envigo`) # 0 significant OTUs

#Combine stats for each individual source for OTUs that varied from day -1 to 1 with adjusted p.values < 0.05----
#Read in individual tables:
for (s in mouse_sources){
  name <- paste("w_sig_otu_", s, sep = "") 
  assign(name, read_tsv(file = paste0("data/process/otu_stats_dn1to0_", s, ".tsv")))
}
#Filter tables so that only significant taxa are shown and make a new column to note source
w_sig_otu_Schloss <- w_sig_otu_Schloss %>% 
filter(p.value.adj < 0.05) %>% 
  mutate(source = "Schloss")
w_sig_otu_Young <- w_sig_otu_Young %>% 
  filter(p.value.adj < 0.05) %>% 
  mutate(source = "Young")
w_sig_otu_Jackson <- w_sig_otu_Jackson %>% 
  filter(p.value.adj < 0.05) %>% 
  mutate(source = "Jackson")
`w_sig_otu_Charles River` <- `w_sig_otu_Charles River` %>% 
  filter(p.value.adj < 0.05) %>% 
  mutate(source = "Charles River")
w_sig_otu_Taconic <- w_sig_otu_Taconic %>% 
  filter(p.value.adj < 0.05) %>% 
  mutate(source = "Taconic")
w_sig_otu_Envigo <- w_sig_otu_Envigo %>% 
  filter(p.value.adj < 0.05) %>% 
  mutate(source = "Envigo")
#Combine tables for all sources of mice
otu_stats_sources_combined <- rbind(w_sig_otu_Schloss, w_sig_otu_Young, w_sig_otu_Jackson, 
                                    `w_sig_otu_Charles River`, w_sig_otu_Taconic, w_sig_otu_Envigo) %>% 
  write_tsv(path = paste0("data/process/otu_dn1to0_stats_all_sources.tsv")) #Save combined dataframe as a .tsv

#OTUs with relative abundances significantly altered by clindamycin treatment that are shared across sources of mice----
#Shared OTUs across all sources of mice:
shared_all_sources <- intersect_all(`sig_otu_Schloss`, `sig_otu_Young`, `sig_otu_Jackson`, `sig_otu_Charles River`, `sig_otu_Taconic`, `sig_otu_Envigo`)
summary(shared_all_sources) # 0 OTUs 

#Shared OTUs across Schloss and Young lab mice:
shared_Schloss_Young <- intersect_all(`sig_otu_Schloss`, `sig_otu_Young`)
summary(shared_Schloss_Young) #2 OTUs 
#"Lactobacillus (OTU 6)", "Lactobacillus (OTU 18)"

#Shared OTUs across Schloss, Young and Charles River mice:
shared_Schloss_Young_CR <- intersect_all(`sig_otu_Schloss`, `sig_otu_Young`, `sig_otu_Charles River`)
summary(shared_Schloss_Young_CR) #0 OTUs 

#Shared OTUs across Jackson, Taconic and Envigo mice:
shared_JAX_Tac_Env <- intersect_all(`sig_otu_Jackson`, `sig_otu_Taconic`, `sig_otu_Envigo`)
summary(shared_JAX_Tac_Env) #0 OTUs

#Shared OTUs across mice purchased from 4 vendors:
shared_4_vendors <- intersect_all(`sig_otu_Jackson`, `sig_otu_Charles River`, `sig_otu_Taconic`, `sig_otu_Envigo`)
summary(shared_4_vendors) #0 OTUss

#Shared otus across Jackson and Charles River mice:
shared_JAX_CR <- intersect_all(`sig_otu_Jackson`, `sig_otu_Charles River`)
summary(shared_JAX_CR) #0 OTUs

#Shared otus across Taconic and Envigo mice:
shared_Tac_Env <- intersect_all(`sig_otu_Taconic`, `sig_otu_Envigo`)
summary(shared_Tac_Env) #0 OTUs

#Function to plot OTUS of interest that overlap with top 20 OTUS in 3 logistic regression models 
#and/or were significantly different across sources of mice at day -1, 0, or 1
#Function to plot specific Otus over time
otu_over_time <- function(otu_plot){
  otu_mean <- agg_otu_data %>% 
    filter(otu == otu_plot) %>% 
    group_by(vendor, day) %>% 
    summarize(mean=(mean(agg_rel_abund + 1/10874))) %>% 
    ungroup
  otu_mice <-  agg_otu_data %>% 
    filter(otu == otu_plot) %>% 
    mutate(agg_rel_abund = agg_rel_abund + 1/10874) %>%
    select(vendor, day, agg_rel_abund, otu)
  otu_time <- ggplot(NULL)+
    geom_point(otu_mice, mapping = aes(x=day, y=agg_rel_abund, color=vendor, alpha = .2), size  = .5, show.legend = FALSE, position = position_dodge(width = 0.6))+
    geom_line(otu_mean, mapping = aes(x=day, y=mean, color=vendor), size = 1)+
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_vendors,
                        labels=color_vendors)+
    geom_hline(yintercept=1/5437, color="gray")+
    labs(title=otu_plot,
         x="Day",
         y="Relative abundance (%)") +
    scale_x_continuous(breaks = c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
                       limits = c(-1.5, 9.5)) +
    scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
    theme_classic()+
    theme(plot.title=element_text(hjust=0.5, face = "italic"))+
    theme(legend.title=element_blank())+
    theme(text = element_text(size = 16))  # Change font size for entire plot
  save_plot(filename = paste0("results/figures/", otu_plot,"_time.png"), otu_time, base_aspect_ratio = 2)
}

#Taxa that differ across mouse sources at multiple timepoints and are important features in 2/3 logistic regresssion models 
otu_over_time("Bacteroides (OTU 2)")
otu_over_time("Enterobacteriaceae (OTU 1)")
otu_over_time("Enterococcus (OTU 23)")

#Taxa that are associated with C. difficile colonization (important to 2/3 logistic regresion models)
otu_over_time("Ruminococcaceae (OTU 520)")
otu_over_time("Erysipelotrichaceae (OTU 234)")
otu_over_time("Porphyromonadaceae (OTU 7)")
otu_over_time("Lactobacillus (OTU 18)")

#Taxa that are different across vendors and are in the top 20 features for the classification model at that particular timepoint
#Day 0 taxa (OTU 1 & 2, which were already plotted above since they appear in multiple models)
otu_over_time("Proteus (OTU 16)")
#Day -1 taxa
otu_over_time("Erysipelotrichaceae (OTU 189)")
otu_over_time("Betaproteobacteria (OTU 58)")
otu_over_time("Burkholderiales (OTU 34)")
otu_over_time("Coriobacteriaceae (OTU 293)")
#Day 1 taxa
otu_over_time("Turicibacter (OTU 14)")
otu_over_time("Bifidobacterium (OTU 46)")
otu_over_time("Bacteroides (OTU 3)")
otu_over_time("Lactobacillus (OTU 31)")
otu_over_time("Bacteroides (OTU 32)")
otu_over_time("Anaerofustis (OTU 475)")

#Additional taxa important to Day 0 logistic regression models 
# but these taxa were not sig. at that timepoint and/or didn't overlap with important taxa from the other models 
otu_over_time("Escherichia/Shigella (OTU 610)") #Error so made without using plot function below
otu_over_time("Lachnospiraceae (OTU 56)")
otu_over_time("Porphyromonadaceae (OTU 54)")
otu_over_time("Lachnospiraceae (OTU 38)")
otu_over_time("Porphyromonadaceae (OTU 22)")
otu_over_time("Lachnospiraceae (OTU 33)")
otu_over_time("Ruminococcaceae (OTU 60)")
otu_over_time("Alishewanella (OTU 776)")
otu_over_time("Lachnospiraceae (OTU 9)")
otu_over_time("Eisenbergiella (OTU 164)")
otu_over_time("Clostridium (OTU 226)")
otu_over_time("Clostridium (OTU 99)")
otu_over_time("Lactobacillus (OTU 834)")

#Need to do separate plot of Escherichia Shigella because the slashes mess up the file name
Escherichia_shigella_over_time <- function(otu_plot){
  otu_mean <- agg_otu_data %>% 
    filter(otu == "Escherichia/Shigella (OTU 610)") %>% 
    group_by(vendor, day) %>% 
    summarize(mean=(mean(agg_rel_abund + 1/10874))) %>% 
    ungroup
  otu_mice <-  agg_otu_data %>% 
    filter(otu == "Escherichia/Shigella (OTU 610)") %>% 
    mutate(agg_rel_abund = agg_rel_abund + 1/10874) %>%
    select(vendor, day, agg_rel_abund, otu)
  otu_time <- ggplot(NULL)+
    geom_point(otu_mice, mapping = aes(x=day, y=agg_rel_abund, color=vendor, alpha = .2), size  = .5, show.legend = FALSE, position = position_dodge(width = 0.6))+
    geom_line(otu_mean, mapping = aes(x=day, y=mean, color=vendor), size = 1)+
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_vendors,
                        labels=color_vendors)+
    geom_hline(yintercept=1/5437, color="gray")+
    labs(title=otu_plot,
         x="Day",
         y="Relative abundance (%)") +
    scale_x_continuous(breaks = c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
                       limits = c(-1.5, 9.5)) +
    scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
    theme_classic()+
    theme(plot.title=element_text(hjust=0.5, face = "italic"))+
    theme(legend.title=element_blank())+
    theme(text = element_text(size = 16))  # Change font size for entire plot
  save_plot(filename = paste0("results/figures/Escherichia_time.png"), otu_time, base_aspect_ratio = 2)
}
Escherichia_shigella_over_time("Escherichia/Shigella (OTU 610)")

#Table combining stats testing for differences across sources of mice for all timepoints
otu_sig_dn1 <- tibble(`sig_otu_day-1`) %>% 
  mutate(day = -1) %>% 
  rename(otu = `sig_otu_day-1`)
otu_sig_d0 <- tibble(`sig_otu_day0`) %>% 
  mutate(day = 0) %>% 
  rename(otu = `sig_otu_day0`)
otu_sig_d1 <- tibble(`sig_otu_day1`) %>% 
  mutate(day = 1) %>% 
  rename(otu = `sig_otu_day1`)
otu_sig_d2 <- tibble(`sig_otu_day2`) %>% 
  mutate(day = 2) %>% 
  rename(otu = `sig_otu_day2`)
otu_sig_d3 <- tibble(`sig_otu_day3`) %>% 
  mutate(day = 3) %>% 
  rename(otu = `sig_otu_day3`)
otu_sig_d4 <- tibble(`sig_otu_day4`) %>% 
  mutate(day = 4) %>% 
  rename(otu = `sig_otu_day4`)
otu_sig_d5 <- tibble(`sig_otu_day5`) %>% 
  mutate(day = 5) %>% 
  rename(otu = `sig_otu_day5`)
otu_sig_d6 <- tibble(`sig_otu_day6`) %>% 
  mutate(day = 6) %>% 
  rename(otu = `sig_otu_day6`)
otu_sig_d7 <- tibble(`sig_otu_day7`) %>% 
  mutate(day = 7) %>% 
  rename(otu = `sig_otu_day7`)
otu_sig_d8 <- tibble(`sig_otu_day8`) %>% 
  mutate(day = 8) %>% 
  rename(otu = `sig_otu_day8`)
otu_sig_d9 <- tibble(`sig_otu_day9`) %>% 
  mutate(day = 9) %>% 
  rename(otu = `sig_otu_day9`)
sig_otus_across_days <- rbind(otu_sig_dn1, otu_sig_d0, otu_sig_d1, otu_sig_d2, otu_sig_d3, otu_sig_d4, otu_sig_d5, otu_sig_d6, otu_sig_d7, otu_sig_d8, otu_sig_d9)



