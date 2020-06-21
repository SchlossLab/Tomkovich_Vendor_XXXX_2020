source("code/functions.R")

#Comparing taxa identified via statistical analysis (differences across sources and changes due to clindamycin treatment) to taxa that showed up as important in logistic regression models----
#Make basic venn diagram without labels to fill with taxa overlap data
df.venn <- data.frame(x = c(10, -10),
                      y = c(-10, -10),
                      labels = c('Source', 'Clindamycin'))
ggplot(df.venn, aes(x0 = x, y0 = y, r = 19, fill = labels)) +
  geom_circle(alpha = .3, size = 1, colour = 'grey') +
  coord_fixed() +
  theme_void()

#Families identified in logistic regression classification models (20 families with the highest ranking for each model)
intrep_families <- read_tsv("data/process/combined_top20_families_all_models.tsv")
interp_families_dn1 <- intrep_families %>% filter(model_input_day == -1) %>% pull(family)
interp_families_d0 <- intrep_families %>% filter(model_input_day == 0) %>% pull(family)
interp_families_d1 <- intrep_families %>% filter(model_input_day == 1) %>% pull(family)
interp_families_combined <- c(interp_families_dn1, interp_families_d0, interp_families_d1)

#Families that vary across sources of mice on day -1, 0, or 1.
sig_family_source <- read_tsv("data/process/family_stats_dn1to1_combined.tsv")
`sig_family_day-1` <- sig_family_source %>% filter(day == -1) %>% pull(family)
`sig_family_day0` <- sig_family_source %>% filter(day == 0) %>% pull(family)
`sig_family_day1` <- sig_family_source %>% filter(day == 1) %>% pull(family)
#Families that consistently vary across sources of mice on day-1, 0, or 1
shared_sig_families_Dn1toD1 <- intersect_all(`sig_family_day-1`, sig_family_day0, sig_family_day1)

#Families that were impacted by clindamycin treatment
`sig_family_pairs` <- read_tsv("data/process/family_stats_dn1to0.tsv") %>% 
  filter(p.value.adj < 0.05) %>% pull(family)

# Day -1 overlapping families----
# Overlap between families that varied across sources of mice on day -1 and top 20 taxa from logistic regression model based on day -1 community:
dayn1_and_interp_dn1_f <- intersect_all(`sig_family_day-1`, `interp_families_dn1`)
# Overlap between families that were altered by clindamycin treatment and top 20 taxa from logistic regression model based on day -1 community:
paired_and_interp_dn1_f <- intersect_all(`sig_family_pairs`, `interp_families_dn1`)

# Day 0 overlapping families----
# Overlap between families that varied across sources of mice on day 0 and top 20 taxa from logistic regression model based on day 0 community:
day0_and_interp_d0_f <- intersect_all(`sig_family_day0`, `interp_families_d0`)
# Overlap between families that were altered by clindamycin treatment and top 20 taxa from logistic regression model based on day 0 community:
paired_and_interp_d0_f <- intersect_all(`sig_family_pairs`, `interp_families_d0`)

# Day 1 overlapping families----
# Overlap between families that varied across sources of mice on day 1 and top 20 taxa from logistic regression model based on day 1 community:
day1_and_interp_d1_f <- intersect_all(`sig_family_day1`, `interp_families_d1`)
# Overlap between families that were altered by clindamycin treatment and top 20 taxa from logistic regression model based on day 1 community:
paired_and_interp_d1_f <- intersect_all(`sig_family_pairs`, `interp_families_d1`)

#Combined overlapping families based on day -1, 0, and 1 comparisons----
#Combined overlapping families that varied by source and were in the top 20 taxa from at least one logistic regression model & remove any duplicates
dayn1to1_and_interp_combined_f <- unique(c(`dayn1_and_interp_dn1_f`, `day0_and_interp_d0_f`, day1_and_interp_d1_f))
#See which source families show up as significant over multiple days
source_families_w_duplicates <- c(`dayn1_and_interp_dn1_f`, `day0_and_interp_d0_f`, day1_and_interp_d1_f)
key_source_families <- sort(unique(source_families_w_duplicates[duplicated(source_families_w_duplicates)]))
#Key source families: "Bacteroidaceae" -seen across 3 days, "Deferribacteraceae", "Enterococcaceae", "Lachnospiraceae"
#Combined overlapping families that were altered by clindamycin treatment and were in the top 20 taxa from at least one logistic regression model
paired_and_interp_combined_f <- unique(c(`paired_and_interp_dn1_f`,`paired_and_interp_d0_f`, paired_and_interp_d1_f))
#See which clindamycin-associated families show up as significant over multiple days
clind_familes_w_duplicates <- c(`paired_and_interp_dn1_f`,`paired_and_interp_d0_f`, paired_and_interp_d1_f)
key_clind_families <- sort(unique(clind_familes_w_duplicates[duplicated(clind_familes_w_duplicates)]))
#Key clindamycin families: "Bifidobacteriaceae"  "Coriobacteriaceae - seen across 3 days, "Enterococcaceae",    
#"Lachnospiraceae", "Ruminococcaceae", "Ruminococcaceae", "Verrucomicrobiaceae"
#List of key families from source and clindamycin varying families:
key_families <- unique(c(key_source_families, key_clind_families))

#Function to plot venn diagrams at the family level
#Arguments:
# source_comp: overlap between variation by source and taxa that were important in classification model
# clind_comp: overlap between variation by source and taxa that were important in classification model
# title: title for plot
venn_families <- function(source_comp, clind_comp, title){
  #Vector of numbers that represent the comparisons between source_comp and clind_comp families
  families <- c(length(setdiff(source_comp, clind_comp)), 
                length(intersect(source_comp, clind_comp)), 
                length(setdiff(clind_comp, source_comp)))  

#Make data frame to annotate number of unique and overlapping taxa onto venn diagram plot
  df_venn_families <- as.data.frame(families) %>%
  mutate(x = c(-15, 0, 15),
         y = c(3, 3, 3))
#List of taxa that are unique to each group or overlap to add to the plot
  source_unique_taxa <- sort(setdiff(source_comp, clind_comp))
# intersection
  clind_unique_taxa <- sort(setdiff(clind_comp, source_comp))  
# intersection
  overlap_taxa <- sort(intersect(source_comp, clind_comp))

#Identify families that overlap with the key_families that show up across multiple comparisons
  source_df <-  tibble(
    families = source_unique_taxa,
    overlap = source_unique_taxa %in% key_families) %>% 
    arrange(families)
  clind_df <-  tibble(
    families = clind_unique_taxa,
    overlap = clind_unique_taxa %in% key_families) %>% 
    arrange(families)
  overlap_df <- tibble(
    families = overlap_taxa,
    overlap = overlap_taxa %in% key_families) %>% 
    arrange(families)
  
  #Annotations for each set of taxa that overlap with key_families or don't overlap for each part of the venn diagram
  source_overlap <- paste(source_df %>% filter(overlap == TRUE) %>% pull(families), collapse="\n")
  source_unique <- paste(source_df %>% filter(overlap == FALSE) %>% pull(families), collapse="\n")
  clind_overlap <- paste(clind_df %>% filter(overlap == TRUE) %>% pull(families), collapse="\n")
  clind_unique <- paste(clind_df %>% filter(overlap == FALSE) %>% pull(families), collapse="\n")
  intersect_overlap <- paste(overlap_df %>% filter(overlap == TRUE) %>% pull(families), collapse="\n")
  intersect_unique <- paste(overlap_df %>% filter(overlap == FALSE) %>% pull(families), collapse="\n")
  
  families_venn_plot <- ggplot(df.venn) +
    geom_circle(aes(x0 = x, y0 = y, r = 19, fill = labels), alpha = .3, size = 1, colour = 'grey') +
    coord_fixed() +
    theme_void() +
    theme(legend.position = 'none') +
    scale_fill_manual(values = c('cornflowerblue', 'firebrick',  'gold')) +
    scale_colour_manual(values = c('cornflowerblue', 'firebrick', 'gold'), guide = FALSE) +
    labs(fill = NULL) +
    annotate("text", x = df_venn_families$x, y = df_venn_families$y, label = df_venn_families[,1], size = 5)+
    annotate("text", x = c(-10, 10), y = c(10, 10), label = c("Source", "Clindamycin"), size = 5)+
    geom_text(label = source_overlap, x = -19, y = -6, size = 2.8, aes(fontface="bold.italic"))+
    geom_text(label = source_unique, x = -19, y = -15, size = 2.8, aes(fontface="italic"))+
    geom_text(label = clind_overlap, x = 19, y = -6, size = 2.8, aes(fontface="bold.italic"))+
    geom_text(label = clind_unique, x = 19, y = -15, size = 2.8, aes(fontface="italic"))+
    geom_text(label = intersect_overlap, x = 0, y = -6, size = 2.8, aes(fontface="bold.italic"))+
    geom_text(label = intersect_unique, x = 0, y = -15, size = 2.8, aes(fontface="italic"))+
    annotate("text", x = 0, y = 14, label = title, size = 5)
}

# Venn diagram of Day -1 overlapping families----
dn1_families_venn_plot <- venn_families(dayn1_and_interp_dn1_f, paired_and_interp_dn1_f, "Day -1 model key family comparisons")
save_plot("results/figures/venn_dn1_families.png", dn1_families_venn_plot, base_aspect_ratio = 1.8)

# Venn diagram of Day 0 overlapping families----
d0_families_venn_plot <- venn_families(day0_and_interp_d0_f, paired_and_interp_d0_f, "Day 0 model key family comparisons")
save_plot("results/figures/venn_d0_families.png", d0_families_venn_plot, base_aspect_ratio = 1.8)

# Venn diagram of Day 1 overlapping families----
d1_families_venn_plot <- venn_families(day1_and_interp_d1_f, paired_and_interp_d1_f, "Day 1 model key family comparisons")
save_plot("results/figures/venn_d1_families.png", d1_families_venn_plot, base_aspect_ratio = 1.8)

# Venn diagram of overlapping families for all days combined----
combined_families_venn_plot <- venn_families(dayn1to1_and_interp_combined_f, paired_and_interp_combined_f, "Key family comparisons for day -1, 0, and 1 models")
save_plot("results/figures/venn_overall_families.png", combined_families_venn_plot, base_aspect_ratio = 1.8)

#OTUs identified in logistic regression classification models (20 OTUs with the highest ranking for each model)
interp_otus <- read_tsv("data/process/combined_top20_otus_all_models.tsv")
interp_otus_dn1 <- interp_otus  %>% filter(model_input_day == -1) %>% pull(OTU)
interp_otus_d0 <- interp_otus  %>% filter(model_input_day == 0) %>% pull(OTU)
interp_otus_d1 <- interp_otus  %>% filter(model_input_day == 1) %>% pull(OTU)
interp_combined <- c(interp_otus_dn1, interp_otus_d0, interp_otus_d1)

#OTUs that vary across sources of mice on day -1, 0, or 1
sig_otu_source <- read_tsv("data/process/otu_stats_dn1to1_combined.tsv")
`sig_otu_day-1` <- sig_otu_source %>% filter(day == -1) %>% pull(otu)
`sig_otu_day0` <- sig_otu_source %>% filter(day == 0) %>% pull(otu)
`sig_otu_day1` <- sig_otu_source %>% filter(day == 1) %>% pull(otu)
#OTUs that consistently vary across sources of mice on day-1, 0, or 1
shared_sig_otus_Dn1toD1 <- intersect_all(`sig_otu_day-1`, sig_otu_day0, sig_otu_day1)


#OTUs that were impacted by clindamycin treatment
`sig_otu_pairs` <- read_tsv("data/process/otu_stats_dn1to0.tsv") %>% 
  filter(p.value.adj < 0.05) %>% pull(otu)

# Day -1 overlapping OTUs----
# Overlap between OTUs that varied across sources of mice on day -1 and top 20 taxa from logistic regression model based on day -1 community:
dayn1_and_interp_dn1_o <- intersect_all(`sig_otu_day-1`, `interp_otus_dn1`)
# Overlap between OTUs that were altered by clindamycin treatment and top 20 taxa from logistic regression model based on day -1 community:
paired_and_interp_dn1_o <- intersect_all(`sig_otu_pairs`, `interp_otus_dn1`)

# Day 0 overlapping OTUs----
# Overlap between OTUs that varied across sources of mice on day 0 and top 20 taxa from logistic regression model based on day 0 community:
day0_and_interp_d0_o <- intersect_all(`sig_otu_day0`, `interp_otus_d0`)
# Overlap between OTUs that were altered by clindamycin treatment and top 20 taxa from logistic regression model based on day 0 community:
paired_and_interp_d0_o <- intersect_all(`sig_otu_pairs`, `interp_otus_d0`)

# Day 1 overlapping OTUs----
# Overlap between OTUs that varied across sources of mice on day 1 and top 20 taxa from logistic regression model based on day 1 community:
day1_and_interp_d1_o <- intersect_all(`sig_otu_day1`, `interp_otus_d1`)
# Overlap between OTUs that were altered by clindamycin treatment and top 20 taxa from logistic regression model based on day 1 community:
paired_and_interp_d1_o <- intersect_all(`sig_otu_pairs`, `interp_otus_d1`)

#Combined overlapping OTUs based on day -1, 0, and 1 comparisons----
#Combined overlapping OTUs that varied by source and were in the top 20 taxa from at least one logistic regression model & remove any duplicates
dayn1to1_and_interp_combined <- unique(c(`dayn1_and_interp_dn1_o`, `day0_and_interp_d0_o`, `day1_and_interp_d1_o`))
#See which source OTUs show up as significant over multiple days
source_otus_w_duplicates <- c(`dayn1_and_interp_dn1_o`, `day0_and_interp_d0_o`, `day1_and_interp_d1_o`)
key_source_otus <- sort(source_otus_w_duplicates[duplicated(source_otus_w_duplicates)])
#Key source OTUs: "Bacteroides (OTU 2)", "Enterococcus (OTU 23)"
#Combined overlapping OTUs that were altered by clindamycin treatment and were in the top 20 taxa from at least one logistic regression model
paired_and_interp_combined <- unique(c(`paired_and_interp_dn1_o`, `paired_and_interp_d0_o`, `paired_and_interp_d1_o`))
#See which clindamycin-associated OTUs show up as significant over multiple days
clind_otus_w_duplicates <- c(`paired_and_interp_dn1_o`, `paired_and_interp_d0_o`, `paired_and_interp_d1_o`)
key_clind_otus <- sort(clind_otus_w_duplicates[duplicated(clind_otus_w_duplicates)])
#Key clindamycin OTUs: "Enterobacteriaceae (OTU 1)", "Enterococcus (OTU 23)", "Porphyromonadaceae (OTU 7)"
#List of key OTUs from source and clindamycin varying OTUs:
key_otus <- unique(c(key_source_otus, key_clind_otus))

#Function to plot venn diagrams at the OTU level
#Arguments:
# source_comp: overlap between variation by source and taxa that were important in classification model
# clind_comp: overlap between variation by source and taxa that were important in classification model
venn_otus <- function(source_comp, clind_comp, title){
  #Vector of numbers that represent the comparisons between source_comp and clind_comp OTUs
  otus <- c(length(setdiff(source_comp, clind_comp)), 
                length(intersect(source_comp, clind_comp)), 
                length(setdiff(clind_comp, source_comp)))
  #Make data frame to annotate number of unique and overlapping OTUs onto venn diagram plot
  df_venn_otus <- as.data.frame(otus) %>%
    mutate(x = c(-15, 0, 15),
           y = c(5, 5, 5))
  #List of OTUs that are unique to each group or overlap to add to the plot
  source_unique_taxa <- sort(setdiff(source_comp, clind_comp))  
  # intersection
  clind_unique_taxa <- sort(setdiff(clind_comp, source_comp))  
  # intersection
  overlap_taxa <- sort(intersect(source_comp, clind_comp))
  
  #Identify families that overlap with the key_families that show up across multiple comparisons
  source_df <-  tibble(
    otus = source_unique_taxa,
    overlap = source_unique_taxa %in% key_otus) %>% 
    arrange(otus)
  clind_df <-  tibble(
    otus = clind_unique_taxa,
    overlap = clind_unique_taxa %in% key_otus) %>% 
    arrange(otus)
  overlap_df <- tibble(
    otus = overlap_taxa,
    overlap = overlap_taxa %in% key_otus) %>% 
    arrange(otus)
  
  #Annotations for each set of taxa that overlap with key_otus or don't overlap for each part of the venn diagram
  source_overlap <- paste(source_df %>% filter(overlap == TRUE) %>% pull(otus), collapse="\n")
  source_unique <- paste(source_df %>% filter(overlap == FALSE) %>% pull(otus), collapse="\n")
  clind_overlap <- paste(clind_df %>% filter(overlap == TRUE) %>% pull(otus), collapse="\n")
  clind_unique <- paste(clind_df %>% filter(overlap == FALSE) %>% pull(otus), collapse="\n")
  intersect_overlap <- paste(overlap_df %>% filter(overlap == TRUE) %>% pull(otus), collapse="\n")
  intersect_unique <- paste(overlap_df %>% filter(overlap == FALSE) %>% pull(otus), collapse="\n")
  
  otus_venn_plot <- ggplot(df.venn) +
    geom_circle(aes(x0 = x, y0 = y, r = 19, fill = labels), alpha = .3, size = 1, colour = 'grey') +
    coord_fixed() +
    theme_void() +
    theme(legend.position = 'none') +
    scale_fill_manual(values = c('cornflowerblue', 'firebrick',  'gold')) +
    scale_colour_manual(values = c('cornflowerblue', 'firebrick', 'gold'), guide = FALSE) +
    labs(fill = NULL) +
    annotate("text", x = df_venn_otus$x, y = df_venn_otus$y, label = df_venn_otus[,1], size = 5)+
    annotate("text", x = c(-10, 10), y = c(10, 10), label = c("Source", "Clindamycin"), size = 5)+
    geom_text(label = source_overlap, x = -19, y = 0.5, size = 2.8, aes(fontface="bold.italic"))+
    geom_text(label = source_unique, x = -19, y = -13, size = 2.8, aes(fontface="italic"))+
    geom_text(label = clind_overlap, x = 19, y = 0.5, size = 2.8, aes(fontface="bold.italic"))+
    geom_text(label = clind_unique, x = 19, y = -13, size = 2.8, aes(fontface="italic"))+
    geom_text(label = intersect_overlap, x = 0, y = 0.5, size = 2.8, aes(fontface="bold.italic"))+
    geom_text(label = intersect_unique, x = 0, y = -13, size = 2.8, aes(fontface="italic"))+
    annotate("text", x = 0, y = 14, label = title, size = 5)
}

# Venn diagram of Day -1 overlapping OTUs----
dn1_otus_venn_plot <- venn_otus(dayn1_and_interp_dn1_o, paired_and_interp_dn1_o, "Day -1 model key OTU comparisons")
save_plot("results/figures/venn_dn1_otus.png", dn1_otus_venn_plot, base_aspect_ratio = 1.8)

# Venn diagram of Day 0 overlapping OTUs----
d0_otus_venn_plot <- venn_otus(day0_and_interp_d0_o, paired_and_interp_d0_o, "Day 0 model key OTU comparisons")
save_plot("results/figures/venn_d0_otus.png", d0_otus_venn_plot, base_aspect_ratio = 1.8)

# Venn diagram of Day -1 overlapping OTUs----
d1_otus_venn_plot <- venn_otus(day1_and_interp_d1_o, paired_and_interp_d1_o, "Day 1 model key OTU comparisons")
save_plot("results/figures/venn_d1_otus.png", d1_otus_venn_plot, base_aspect_ratio = 1.8)

# Venn diagram of combined overlapping OTUs----
combined_otus_venn_plot <- venn_otus(dayn1to1_and_interp_combined, paired_and_interp_combined, "Key taxa comparisons for day -1, 0, and 1 models")
save_plot("results/figures/venn_overall_otus.png", combined_otus_venn_plot, base_aspect_ratio = 1.8)

#Comparison of combined Venn diagram OTUs to significant OTUs that varied by source on day -1, 0, and 1:
intersect_all(`shared_sig_otus_Dn1toD1`, `dayn1to1_and_interp_combined`)
intersect_all(`shared_sig_otus_Dn1toD1`, `paired_and_interp_combined`)

#Comparison of combined Venn diagram families to significant families that varied by source on day -1, 0, and 1:
intersect_all(`shared_sig_families_Dn1toD1`, `dayn1to1_and_interp_combined_f`)
intersect_all(`shared_sig_families_Dn1toD1`, `paired_and_interp_combined_f`)
