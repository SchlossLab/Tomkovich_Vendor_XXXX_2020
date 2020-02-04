source("code/functions.R")

# Import taxonomy into data frame and clean up taxonomy names
taxonomy <- read_tsv(file="data/process/vendors.taxonomy") %>%
  rename_all(tolower) %>% #remove uppercase from column names
  # Split taxonomic information into separate columns for each taxonomic level  
  mutate(taxonomy=str_replace_all(taxonomy, c("\\(\\d*\\)" = "", #drop digits with parentheses around them
                                              ';$' = "", #removes semi-colon at end of line
                                              'Bacteria_unclassified' = 'Unclassified',
                                              "Clostridium_" = "Clostridium ", #Remove underscores after Clostridium
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

#Make sure weight & C. diff CFU data are included in metadata frame----
source("code/weight_plot.R") #ERROR, not sourcing correctly. Move to functions
source("code/cfu_plot.R") ##ERROR, not sourcing correctly. Move to functions
#Select only relevant columns from weight & cfu dataframes, plus the id column for merging dataframes
percent_baseline_weight_data <- percent_baseline_weight_data %>% 
  ungroup() %>% 
  select(id, baseline_weight, percent_baseline_weight, lowest_percent_baseline_weight)
cfu_data_final <- cfu_data_final %>% 
  ungroup() %>% 
  select(id, cfu, cfu_d3, cfu_d4, cfu_d5, cfu_d6, cfu_d7)
#Join weight data to metadata
metadata <- left_join(metadata, percent_baseline_weight_data, by = "id")
#Join cfu data to metadata
metadata <- left_join(metadata, cfu_data_final, by = "id") 

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

#Kruskal_wallis test for differences across groups at different taxonomic levels with Benjamini-Hochburg correction
#Arguments:
# dataframe=dataframe to analyze
# timepoint = timepoint being analyzed
#taxonomic_level = taxon (family, genus, etc)
kruskal_wallis_groups <- function(dataframe, timepoint, taxonomic_level){
  taxonomic_level <- enquo(taxonomic_level)
  dataframe %>% 
    filter(day == timepoint) %>% 
    group_by(!!taxonomic_level) %>% 
    do(tidy(kruskal.test(agg_rel_abund~factor(vendor), data=.))) %>% ungroup() %>% 
    mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
    arrange(p.value.adj) 
}

#Kruskal_wallis test for family differences across sources of mice with Benjamini-Hochburg correction----
#Analyze the following days of the experiment: -1, 0, 1, 2, 5, 8, 9 
#Creates a family_tests_day_ data frame for each date of interest.
#Also creates a sig_family_dayX list of all the significant familes for each date of interest. Significance based on adjusted p value < 0.05.
dates_of_interest <- c(-1, 0, 1, 2, 5, 8, 9)
for(d in dates_of_interest){
  name <- paste("family_tests_day", d, sep = "") #Way to name the data frames based on the date of interest
  sig_name <- paste("sig_family_day", d, sep ="") 
  assign(sig_name, pull_significant_taxa(assign(name, kruskal_wallis_groups(agg_family_data, d, family)), family))
}

#Kruskal_wallis test for genus differences across sources of mice with Benjamini-Hochburg correction---- 
#Analyze the following days of the experiment: -1, 0, 1, 2, 5, 8, 9 
#Creates a genus_tests_day_ data frame for each date of interest.
#Also creates a sig_genus_dayX list of all the significant genera for each date of interest. Significance based on adjusted p value < 0.05.
for(d in dates_of_interest){
  name <- paste("genus_tests_day", d, sep = "") #Way to name the data frames based on the date of interest
  sig_name <- paste("sig_genus_day", d, sep ="") 
  assign(sig_name, pull_significant_taxa(assign(name, kruskal_wallis_groups(agg_genus_data, d, genus)), genus))
}

#Only 9 significant genera at D2. Check to see if we have a lower amount of samples on D2
D2_samples <- metadata %>% filter(day == 2) 
#49 samples
D2_sequence_data <- agg_genus_data %>% filter(day == 2,
                                              genus == "Clostridium sensu_stricto")
#Only 18 samples with sequence data for D2.

#For significant genera at each relevant timepoint, do pairwise.wilcox.test to determine which sources of mice are significantly different from each other.
pairwise.wilcox_groups <- function(timepoint, genera){
  taxa_data <- agg_genus_data %>% 
    filter(day == timepoint,
           genus == genera)
  tidy(pairwise.wilcox.test(g = taxa_data$vendor, x = taxa_data$agg_rel_abund, p.adjust.method = "BH"))
}

# Do pairwise.wilcox tests with BH correction for D-1 timepoint----
for(g in `sig_genus_day-1`){
  name <- paste("pairwise_wilcox_day-1_", g, sep = "") #Way to name the data frames based on the date of interest
  assign(name, pairwise.wilcox_groups(-1, g))
}
# To do: pairwise.wilcox tests for rest of timepoints. How to output these to a table

#Function to find which significant genera/families are shared across days
intersect_all <- function(a,b,...){
  Reduce(intersect, list(a,b,...))
}

#Shared significant genera across days----
shared_sig_genera <- intersect_all(`sig_genus_day-1`, sig_genus_day0, sig_genus_day1, sig_genus_day2, sig_genus_day5, sig_genus_day8, sig_genus_day9)
#Betaproteobacteria Unclassified, Sutterellaceae Unclassified, Parasutterella, Mucispirillum

#Shared significant genera across days, but excluding day 2----
shared_sig_genera_wo_D2 <- intersect_all(`sig_genus_day-1`, sig_genus_day0, sig_genus_day1, sig_genus_day5, sig_genus_day8, sig_genus_day9)
#Betaproteobacteria Unclassified, Sutterellaceae Unclassified, Parasutterella, Bacteroides, Mucispirillum, Turicibacter, Clostridium XVIII, Proteus

#Shared significant genera across D-1 to D1----
shared_sig_genera_Dn1toD1 <- intersect_all(`sig_genus_day-1`, sig_genus_day0, sig_genus_day1)
#Betaproteobacteria Unclassified, Sutterellaceae Unclassified, Parasutterella, Bacteroides, Mucispirillum, Turicibacter, Clostridium XVIII, Enterococcus, Proteus

#Shared significant familes across days, but excluding day 2----
shared_sig_families <- intersect_all(`sig_family_day-1`, sig_family_day0, sig_family_day1, sig_family_day2, sig_family_day5, sig_family_day8, sig_family_day9)
# Betaproteobacteria Unclassified, Sutterellaceae, Deferribacteraceae

#Shared significant familes across days, but excluding day 2----
shared_sig_families_wo_D2 <- intersect_all(`sig_family_day-1`, sig_family_day0, sig_family_day1, sig_family_day5, sig_family_day8, sig_family_day9)
# Betaproteobacteria Unclassified, Bacteroidaceae, Sutterellaceae, Deferribacteraceae

#Kruskal_wallis test for differences within a single source of mice across time at different taxonomic levels with Benjamini-Hochburg correction
#Arguments:
# dataframe=dataframe to analyze (agg_genus_data or agg_taxa_data filtered by day (only include key days in the experiment -1, 0, 1, 3, 4, 5, 8, 9). Also, need to make sure day is treated as a factor in the input data frame
# mouse_source (aka vendor) = the group of mice being analyzed (Schloss, Young, Charles River, Jackson, Taconic or Envigo)
#taxonomic_level = taxon (family, genus, etc)
kruskal_wallis_time <- function(dataframe, mouse_source, taxonomic_level){
  taxonomic_level <- enquo(taxonomic_level)
  dataframe %>% 
    filter(vendor == mouse_source) %>% 
    group_by(!!taxonomic_level) %>% 
    do(tidy(kruskal.test(agg_rel_abund~factor(day), data=.))) %>% ungroup() %>% 
    mutate(p.value.adj=p.adjust(p.value, method="BH")) %>% 
    arrange(p.value.adj) 
}

#Customize agg_taxa_data dataframes to only select the dates of interest and turn day into a factor
agg_genus_data_limited_days <- agg_genus_data %>% 
  filter(day %in% c(-1, 0, 1, 3, 4, 5, 8, 9)) %>% 
  mutate(day=factor(day, levels=c("-1", "0", "1", "3", "4", "5", "8", "9")))

agg_family_data_limited_days <- agg_family_data %>% 
  filter(day %in% c(-1, 0, 1, 3, 4, 5, 8, 9)) %>% 
  mutate(day=factor(day, levels=c("-1", "0", "1", "3", "4", "5", "8", "9")))

#Kruskal_wallis test for genus differences across time within a single colony source of mice with Benjamini-Hochburg correction---- 
#Analyze the following days of the experiment: -1, 0, 1, 3, 4, 5, 8, 9 
#Creates a genus_tests_source data frame for each date of interest.
#Also creates a sig_genus_source list of all the significant genera for each date of interest. Significance based on adjusted p value < 0.05.
mouse_sources <- levels(metadata$vendor)
for(s in mouse_sources){
  name <- paste("genus_tests_", s, sep = "") #Way to name the data frames based on the date of interest
  sig_name <- paste("sig_genus_", s, sep ="") 
  assign(sig_name, pull_significant_taxa(assign(name, kruskal_wallis_time(agg_genus_data_limited_days, s, genus)), genus))
}

# Number of significant genera for each source of mice:
summary(`sig_genus_Schloss`) # 26 significant genera
summary(`sig_genus_Young`) # 29 significant genera
summary(`sig_genus_Jackson`) # 32 significant genera
summary(`sig_genus_Charles River`) # 20 significant genera
summary(`sig_genus_Taconic`) # 7 significant genera
summary(`sig_genus_Envigo`) # 17 significant genera

#Shared significant genera (identified by comparing across timepoints within a specific source of mice) across groups of mice----
#Shared genera across all sources of mice:
shared_all_sources <- intersect_all(`sig_genus_Schloss`, `sig_genus_Young`, `sig_genus_Jackson`, `sig_genus_Charles River`, `sig_genus_Taconic`, `sig_genus_Envigo`)
summary(shared_all_sources) # 1 genera that significantly change over time our shared between Schloss & Young mice
#Peptostreptococcaceae Unclassified is the 1 genus that is significantly changing over time and shared across all sources of mice. This OTU is probably C. difficile so this makes sense

#Shared genera across Schloss and Young lab mice:
shared_Schloss_Young <- intersect_all(`sig_genus_Schloss`, `sig_genus_Young`)
summary(shared_Schloss_Young) #24 genera that significantly change over time are shared between Schloss & Young mice

#Shared genera across Schloss, Young and Charles River mice:
shared_Schloss_Young_CR <- intersect_all(`sig_genus_Schloss`, `sig_genus_Young`, `sig_genus_Charles River`)
summary(shared_Schloss_Young_CR) #12 genera that significantly change over time are shared between Schloss, Young, and Charles River mice. These 3 sources grouped together (cleared faster, less weight loss) on C. diff CFU & Weight loss plots that combined the 2 experiments
#Includes Enterobacteriaceae Unclassified, Clostridiales Unclassified, Acetatifactor, Lachnospiraceae Unclassified, Oscillibacter, Akkermansia, Ruminococcaceae Unclassified, Firmicutes Unclassified, Clostridium XlVb, Pseudoflavonifractor, Unclassified  

#Shared genera across Jackson, Taconic and Envigo mice:
shared_JAX_Tac_Env <- intersect_all(`sig_genus_Jackson`, `sig_genus_Taconic`, `sig_genus_Envigo`)
summary(shared_JAX_Tac_Env) #3 genera that significantly change over time are shared between Jackson, Taconic, and Envigo mice. These 3 sources grouped together (slower clearance, more weight loss) on C. diff CFU & Weight loss plots that combined the 2 experiments
#Includes Peptostreptococcaceae Unclassified, Enterococcus, Coriobacteriaceae Unclassified

#Shared genera across mice purchased from 4 vendors:
shared_4_vendors <- intersect_all(`sig_genus_Jackson`, `sig_genus_Charles River`, `sig_genus_Taconic`, `sig_genus_Envigo`)
summary(shared_4_vendors) #2 genera Peptostreptococcus and Enterococcus

#Shared genera across Jackson and Charles River mice:
shared_JAX_CR <- intersect_all(`sig_genus_Jackson`, `sig_genus_Charles River`)
summary(shared_JAX_CR) #16 genera

#Shared genera across Taconic and Envigo mice:
shared_Tac_Env <- intersect_all(`sig_genus_Taconic`, `sig_genus_Envigo`)
summary(shared_Tac_Env) #4 genera Peptostreptococcaceae Unclassified, Enterococcus, Coriobacteriaceae Unclassified, Parabacteroides

#Kruskal_wallis test for family differences across time within a single colony source of mice with Benjamini-Hochburg correction---- 
#Analyze the following days of the experiment: -1, 0, 1, 3, 4, 5, 8, 9 
#Creates a family_tests_source data frame for each date of interest.
#Also creates a sig_family_source list of all the significant families for each date of interest. Significance based on adjusted p value < 0.05.
for(s in mouse_sources){
  name <- paste("family_tests_", s, sep = "") #Way to name the data frames based on the date of interest
  sig_name <- paste("sig_family_", s, sep ="") 
  assign(sig_name, pull_significant_taxa(assign(name, kruskal_wallis_time(agg_family_data_limited_days, s, family)), family))
}

#Shared significant families (identified by comparing across timepoints within a specific source of mice) across groups of mice----
#Shared families across all sources of mice:
shared_all_sources_families <- intersect_all(`sig_family_Schloss`, `sig_family_Young`, `sig_family_Jackson`, `sig_family_Charles River`, `sig_family_Taconic`, `sig_family_Envigo`)
summary(shared_all_sources_families) # 1 family that significantly changed over time our shared between Schloss & Young mice
#Enterobacteriaceae is the family that is significantly changing over time and shared across all sources of mice. This OTU is probably C. difficile so this makes sense

#Shared families across Schloss and Young lab mice:
shared_Schloss_Young_families <- intersect_all(`sig_family_Schloss`, `sig_family_Young`)
summary(shared_Schloss_Young_families) #15 families that significantly change over time are shared between Schloss & Young mice

#Shared families across Schloss, Young and Charles River mice:
shared_Schloss_Young_CR_families <- intersect_all(`sig_family_Schloss`, `sig_family_Young`, `sig_family_Charles River`)
summary(shared_Schloss_Young_CR_families) #8 families that significantly change over time are shared between Schloss, Young, and Charles River mice. These 3 sources grouped together (cleared faster, less weight loss) on C. diff CFU & Weight loss plots that combined the 2 experiments
#Includes Enterobacteriaceae, Clostridiales Unclassified, Lachnospiraceae, Coriobacteriaceae, Verrucomicrobiaceae, Ruminococcaceae, Firmicutes Unclassified, Unclassified


#Shared families across Jackson, Taconic and Envigo mice:
shared_JAX_Tac_Env_families <- intersect_all(`sig_family_Jackson`, `sig_family_Taconic`, `sig_family_Envigo`)
summary(shared_JAX_Tac_Env_families) #3 families that significantly change over time are shared between Jackson, Taconic, and Envigo mice. These 3 sources grouped together (slower clearance, more weight loss) on C. diff CFU & Weight loss plots that combined the 2 experiments
#Includes Peptostreptococcaceae, Enterococcaceae, Enterobacteriaceae

#Shared families across mice purchased from 4 vendors:
shared_4_vendors_families <- intersect_all(`sig_family_Jackson`, `sig_family_Charles River`, `sig_family_Taconic`, `sig_family_Envigo`)
summary(shared_4_vendors_families) #2 families 

#Shared families across Jackson and Charles River mice:
shared_JAX_CR_families <- intersect_all(`sig_family_Jackson`, `sig_family_Charles River`)
summary(shared_JAX_CR_families) #11 families

#Shared families across Taconic and Envigo mice:
shared_Tac_Env_families <- intersect_all(`sig_family_Taconic`, `sig_family_Envigo`)
summary(shared_Tac_Env_families) #4 families 


# Plots----

#Function to plot all significant genus relative abundances across vendors at a specific timepoint----
#Arguments:
# list_of_genera = list of genera to plot. Example: sig_genus_day2
# timepoint = timepoint to be analyzed
plot_genera_timepoint <- function(list_of_genera, timepoint){
  plot_genera <- agg_genus_data %>% 
    filter(genus %in% list_of_genera) %>% 
    filter(day == timepoint) %>% 
    mutate(genus=factor(genus, list_of_genera)) %>% 
    mutate(agg_rel_abund = agg_rel_abund + 1/10874) %>% # 10,874 is 2 times the subsampling parameter of 5437
    ggplot(aes(x= reorder(genus, agg_rel_abund), y=agg_rel_abund, color=vendor))+
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_vendors,
                        labels=color_vendors)+
    geom_hline(yintercept=1/5437, color="gray")+
    geom_boxplot(outlier.shape = NA, size = 1.2)+
    geom_jitter(shape=19, size=2, alpha=0.6, position=position_jitterdodge(dodge.width=0.7, jitter.width=0.2)) +
    labs(title=NULL, 
         x=NULL,
         y="Relative abundance (%)")+
    scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
    coord_flip()+
    theme_classic()+
    theme(axis.text.y = element_text(face = "italic"))+ #Have the genera show up as italics
    theme(plot.title=element_text(hjust=0.5))+
    theme(legend.position = c(0.85, 0.2)) + #Get rid of legend title & move legend position
    theme(text = element_text(size = 16))  # Change font size for entire plot
  save_plot(filename = paste0("results/figures/genera_assoc_w_", timepoint,".png"), plot_genera, base_height = 11, base_width = 8.5, base_aspect_ratio = 2)
}

#Plot of significant genera at day -1----
plot_genera_timepoint(`sig_genus_day-1`, -1)
#Plot of significant genera at day 0----
plot_genera_timepoint(sig_genus_day0, 0)
#Plot of significant genera at day 1----
plot_genera_timepoint(sig_genus_day1, 1)
#Plot of significant genera at day 2----
plot_genera_timepoint(sig_genus_day2, 2)
#Plot of significant genera at day 5----
plot_genera_timepoint(sig_genus_day5, 5)
#Plot of significant genera at day 8----
plot_genera_timepoint(sig_genus_day8, 8)
#Plot of significant genera at day 9----
plot_genera_timepoint(sig_genus_day9, 9)


#Function to plot 1 significant genus relative abundances across vendors at a specific timepoint----
#Arguments:
# name = name of genus to plot. Example: Enterococcus
# timepoint = timepoint to be analyzed
plot_genus_timepoint <- function(name, timepoint){
  plot_genera <- agg_genus_data %>% 
    filter(genus == name) %>% 
    filter(day == timepoint) %>% 
    mutate(genus=factor(genus, name)) %>% 
    mutate(agg_rel_abund = agg_rel_abund + 1/10874) %>% # 10,874 is 2 times the subsampling parameter of 5437
    ggplot(aes(x= reorder(genus, agg_rel_abund), y=agg_rel_abund, color=vendor))+
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
    theme_classic()+
    theme(axis.text.y = element_text(face = "italic"))+ #Have the genera show up as italics
    theme(plot.title=element_text(hjust=0.5))+
    theme(legend.position = "bottom") + #Get rid of legend title & move legend position
    theme(text = element_text(size = 16))  # Change font size for entire plot
  save_plot(filename = paste0("exploratory/notebook/day", timepoint, "/", name, "_at_day", timepoint,".png"), plot_genera, base_height = 11, base_width = 8.5, base_aspect_ratio = 2)
}

#Plot all significant genera at Day-1 timepoint separately
for(g in `sig_genus_day-1`){
  plot_genus_timepoint(g, -1)  
}
#Plot all significant genera at Day0 timepoint separately
for(g in `sig_genus_day0`){
  plot_genus_timepoint(g, 0)  
}
#Plot all significant genera at Day1 timepoint separately
for(g in `sig_genus_day1`){
  plot_genus_timepoint(g, 1)  
}
#Plot all significant genera at Day2 timepoint separately
for(g in `sig_genus_day2`){
  plot_genus_timepoint(g, 2)  
}
# This loop isn't working because of the slash in Escherichia/Shigella. Temporary fix for now:
sig_genus_day2 <- str_replace(sig_genus_day2, "Escherichia/Shigella", "Escherichia-Shigella")
for(g in `sig_genus_day2`){
  plot_genus_timepoint(g, 2)  
}
#To Do: This generates an empty graph for Escherichia/Shigella, so repeat plot_genus_timepoint but make specific to Escherichia/Shigella

#Plot all significant genera at Day5 timepoint separately
for(g in `sig_genus_day5`){
  plot_genus_timepoint(g, 5)  
}
#Plot all significant genera at Day8 timepoint separately
for(g in `sig_genus_day8`){
  plot_genus_timepoint(g, 8)  
}
#Plot all significant genera at Day9 timepoint separately
for(g in `sig_genus_day9`){
  plot_genus_timepoint(g, 9)  
}

#Function to plot family relative abundances across vendors at a specific timepoint
#Arguments:
# list_of_families = list of families to plot. Example: sig_family_day2
# timepoint = timepoint to be analyzed
plot_families_timepoint <- function(list_of_families, timepoint){
  plot_families <- agg_family_data %>% 
    filter(family %in% list_of_families) %>% 
    filter(day == timepoint) %>% 
    mutate(family=factor(family, list_of_families)) %>% 
    mutate(agg_rel_abund = agg_rel_abund + 1/10874) %>% # 10,874 is 2 times the subsampling parameter of 5437
    ggplot(aes(x= reorder(family, agg_rel_abund), y=agg_rel_abund, color=vendor))+
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_vendors,
                        labels=color_vendors)+
    geom_hline(yintercept=1/5437, color="gray")+
    geom_boxplot(outlier.shape = NA, size = 1.2)+
    geom_jitter(shape=19, size=2, alpha=0.6, position=position_jitterdodge(dodge.width=0.7, jitter.width=0.2)) +
    labs(title=NULL, 
         x=NULL,
         y="Relative abundance (%)")+
    scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
    coord_flip()+
    theme_classic()+
    theme(axis.text.y = element_text(face = "italic"))+ #Have the genera show up as italics
    theme(plot.title=element_text(hjust=0.5))+
    theme(legend.position = c(0.85, 0.2)) + #Get rid of legend title & move legend position
    theme(text = element_text(size = 16))  # Change font size for entire plot
  save_plot(filename = paste0("results/figures/families_assoc_w_", timepoint,".png"), plot_families, base_height = 11, base_width = 8.5, base_aspect_ratio = 2)
}

#Plot of significant families at day -1----
plot_families_timepoint(`sig_family_day-1`, -1)
#Plot of significant families at day 0----
plot_families_timepoint(sig_family_day0, 0)
#Plot of significant families at day 1----
plot_families_timepoint(sig_family_day1, 1)
#Plot of significant families at day 2----
plot_families_timepoint(sig_family_day2, 2)
#Plot of significant families at day 5----
plot_families_timepoint(sig_family_day5, 5)
#Plot of significant families at day 8----
plot_families_timepoint(sig_family_day8, 8)
#Plot of significant families at day 9----
plot_families_timepoint(sig_family_day9, 9)

#Function to plot specific genera over time
genera_over_time <- function(genus_plot){
#  genus_plot <- enquo(genus_plot)
  genus_mean <- agg_genus_data %>% 
    filter(genus == genus_plot) %>% 
    group_by(vendor, day) %>% 
    summarize(mean=(mean(agg_rel_abund + 1/10874))) %>% 
    ungroup
  genus_mice <-  agg_genus_data %>% 
    filter(genus == genus_plot) %>% 
    mutate(agg_rel_abund = agg_rel_abund + 1/10874) %>%
    select(vendor, day, agg_rel_abund, genus)
  genus_time <- ggplot(NULL)+
    geom_point(genus_mice, mapping = aes(x=day, y=agg_rel_abund, color=vendor, alpha = .2), show.legend = FALSE, size = 2.5)+
    geom_line(genus_mean, mapping = aes(x=day, y=mean, color=vendor), size = 1)+
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_vendors,
                        labels=color_vendors)+
    geom_hline(yintercept=1/5437, color="gray")+
    labs(title=genus_plot,
         x="Day",
         y="Relative abundance (%)") +
    scale_x_continuous(breaks = c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
                       limits = c(-1.5, 9.5)) +
    scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
    theme_classic()+
    theme(plot.title=element_text(hjust=0.5, face = "italic"))+
    theme(legend.title=element_blank())+
    theme(text = element_text(size = 16))  # Change font size for entire plot
  save_plot(filename = paste0("results/figures/", genus_plot,"_time.png"), genus_time, base_aspect_ratio = 2)
}

#Plot of C. difficile Otu over time across all timepoints
genera_over_time("Peptostreptococcaceae Unclassified") #Otu0020 corresponds to C. difficile? 

#Plots of shared significant genera across D-1 to D1 for all timepoints----
shared_sig_genera_Dn1toD1 #= Betaproteobacteria Unclassified, Sutterellaceae Unclassified, Parasutterella, Bacteroides, Mucispirillum, Turicibacter, Clostridium XVIII, Enterococcus, Proteus
genera_over_time("Enterococcus")
genera_over_time("Betaproteobacteria Unclassified")
genera_over_time("Sutterellaceae Unclassified")
genera_over_time("Parasutterella")
genera_over_time("Bacteroides")
genera_over_time("Mucispirillum")
genera_over_time("Turicibacter")
genera_over_time("Clostridium XVIII")
genera_over_time("Proteus")

#Function to plot specific familes over time
family_over_time <- function(family_plot){
  #  family_plot <- enquo(family_plot)
  family_mean <- agg_family_data %>% 
    filter(family == family_plot) %>% 
    group_by(vendor, day) %>% 
    summarize(mean=(mean(agg_rel_abund + 1/10874))) %>% 
    ungroup
  family_mice <-  agg_family_data %>% 
    filter(family == family_plot) %>% 
    mutate(agg_rel_abund = agg_rel_abund + 1/10874) %>%
    select(vendor, day, agg_rel_abund, family)
  family_time <- ggplot(NULL)+
    geom_point(family_mice, mapping = aes(x=day, y=agg_rel_abund, color=vendor, alpha = .2), show.legend = FALSE, size = 2.5)+
    geom_line(family_mean, mapping = aes(x=day, y=mean, color=vendor), size = 1)+
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_vendors,
                        labels=color_vendors)+
    geom_hline(yintercept=1/5437, color="gray")+
    labs(title=family_plot,
         x="Day",
         y="Relative abundance (%)") +
    scale_x_continuous(breaks = c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
                       limits = c(-1.5, 9.5)) +
    scale_y_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
    theme_classic()+
    theme(plot.title=element_text(hjust=0.5, face = "italic"))+
    theme(legend.title=element_blank())+
    theme(text = element_text(size = 16))  # Change font size for entire plot
  save_plot(filename = paste0("results/figures/", family_plot,"_time.png"), family_time, base_aspect_ratio = 2)
}
#Plot Clostridiaceae_1 family over time. Could include SFB, but unlikely because it is primarily only in Jackson mice while SFB has been associated primarily with Taconic mice----
family_over_time("Clostridiaceae_1")

# Exploratory correlation analyses----

#Does C. difficile cfu correlate with max. amount of weight lost----
metadata_day1 <- metadata %>% filter(day == 1)
C.diff_weight_corr_1 <- cor.test(metadata_day1$cfu, metadata_day1$lowest_percent_baseline_weight, method = "spearman")
#rho = 0.181, p-value = 0.2622
metadata_day5 <- metadata %>% filter(day == 5)
C.diff_weight_corr_5 <- cor.test(metadata_day5$cfu, metadata_day5$lowest_percent_baseline_weight, method = "spearman")
#rho = -0.509, p-value = 0.0005671

# Plot of cfu at day 5 versus max. amount of weight lost----
cfu_weight_corr_day5_lowest_weight <- ggplot(metadata_day5)+
  geom_point(mapping = aes(x=cfu, y=lowest_percent_baseline_weight, color=vendor, alpha = .2), show.legend = FALSE, size = 2.5)+
  scale_colour_manual(name=NULL,
                      values=color_scheme,
                      breaks=color_vendors,
                      labels=color_vendors)+
  geom_hline(yintercept=1/5437, color="gray")+
  scale_x_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
  labs(title= "Day 5 CFU versus Lowest Percent Baseline Weight",
       x="C. difficile CFU",
       y="Lowest Percent Baseline Weight") +
  theme_classic()+
  theme(plot.title=element_text(hjust=0.5, face = "italic"))+
  theme(legend.title=element_blank())
#  theme(text = element_text(size = 16))  # Change font size for entire plot
save_plot("results/figures/cfu_weight_corr_day5_lowest_weight.png", cfu_weight_corr_day5_lowest_weight, base_aspect_ratio = 2)

#Does C. difficile cfu correlate with amount of weight lost on day 5----
C.diff_weight_corr_5 <- cor.test(metadata_day5$cfu, metadata_day5$percent_baseline_weight, method = "spearman")
#rho = 0.405, p-value = 0.006932

#Are there genera that significantly correlate with the amount of weight lost at specific timepoints?----
for (g in sig_genus_day0){
  agg_genus_data %>% filter(day == 0) %>% filter(genus == g)
  name <- paste("Spearman_day0_cor_", g, sep ="")
  assign(name, cor.test(agg_genus_data$agg_rel_abund, agg_genus_data$percent_baseline_weight, method ="spearman"))
}
#These significant genera vary across mouse sources, may not be related to C. difficile infection dynamics.
#Next step find genera that are consistently varying over key timepoints of the 9 day C. difficile infection.
#See which genera significantly vary over time and are shared across the colony sources. See if these correlate with C. difficile relative abundance at a particular timepoint.

for (g in sig_genus_day1){
  agg_genus_data %>% filter(day == 1) %>% filter(genus == g)
  name <- paste("Spearman_day1_cor_", g, sep ="")
  assign(name, cor.test(agg_genus_data$agg_rel_abund, agg_genus_data$percent_baseline_weight, method ="spearman"))
}

#Alternative: pick timepoints and taxa of interest based on above tests on all family/genus identified from Kruskal-Wallis test with Benjamini-Hochburg correction above
#See if the relative abundances of significant taxa correlate with C. difficile relative abundance or amount of weight lost at a specific timepoint----

Enterococcous_data_D0 <- agg_genus_data %>% filter(genus == "Enterococcus" & day == 0)
Enterococcus_corr_D0_lowest_weight <- cor.test(Enterococcous_data$agg_rel_abund, Enterococcous_data$lowest_percent_baseline_weight, method = "spearman")
#Not significant p-value = 0.1535

Enterococcus_D0_corr_C.diff_D5 <- cor.test(Enterococcous_data_D0$agg_rel_abund, Enterococcous_data_D0$cfu_d5, method = "spearman")
#Significant: p-value = 0.00102; rho = 0.52

#Test Enterobacteriaceae correlations
Enterobacteriaceae_data_D0 <- agg_family_data %>% filter(family == "Enterobacteriaceae" & day == 0) 
Enterobacteriaceae_D0_corr_C.diff_D5 <- cor.test(Enterobacteriaceae_data_D0$agg_rel_abund, Enterobacteriaceae_data_D0$cfu_d5, method = "spearman")
# p = 0.0224, rho = -0.38

# Plot of cfu at day 5 versus Enterococcus relative abundance D0----
cfu_d5_EnterococcusD0 <- ggplot(Enterococcous_data_D0)+
  geom_point(mapping = aes(x=cfu_d5, y=agg_rel_abund, color=vendor, alpha = .2), show.legend = FALSE, size = 2.5)+
  scale_colour_manual(name=NULL,
                      values=color_scheme,
                      breaks=color_vendors,
                      labels=color_vendors)+
  geom_hline(yintercept=1/5437, color="gray")+
  scale_x_log10(breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels=c(1e-2, 1e-1, 1, 10, 100))+
  labs(title= "Enterococcus D0 versus C. difficile CFU D5",
       x="C. difficile CFU D5",
       y="Enterococcus D0 relative abundance") +
  theme_classic()+
  theme(plot.title=element_text(hjust=0.5, face = "italic"))+
  theme(legend.title=element_blank())
#  theme(text = element_text(size = 16))  # Change font size for entire plot
save_plot("results/figures/cfu_d5_EnterococcusD0.png", cfu_d5_EnterococcusD0, base_aspect_ratio = 2)
