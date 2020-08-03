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

#OTUs identified in logistic regression classification models (20 OTUs with the highest ranking for each model)
interp_otus <- read_tsv("data/process/combined_top20_otus_all_models.tsv")

#Function to get top 20 OTUs from each model:
get_interp_otus <- function(timepoint){
  otus <- interp_otus  %>% filter(model_input_day == timepoint) %>% 
    separate(OTU, into = c("bactname", "OTUnumber"), sep = "\\ [(]", remove = FALSE) %>% #Add columns to separate bacteria name from OTU number to utilize ggtext so that only bacteria name is italicized
    mutate(otu_name = glue("*{bactname}* ({OTUnumber}")) %>% #Markdown notation so that only bacteria name is italicized
    pull(otu_name)
}
#format_otus <- (paste(otus, collapse = "<br>")) #Format so that there's line breaks separating each OTU
#return(format_otus)

interp_otus_dn1 <- get_interp_otus(-1)
interp_otus_d0 <- get_interp_otus(0)
interp_otus_d1 <- get_interp_otus(1)
interp_combined <- c(interp_otus_dn1, interp_otus_d0, interp_otus_d1)

#OTUs that vary across sources of mice on day -1, 0, or 1
sig_otu_source <- read_tsv("data/process/otu_stats_dn1to1_combined.tsv")
#Function to get OTUs that vary across sources for each timepoint:
get_source_otus <- function(timepoint){
  otus <- sig_otu_source  %>% filter(day == timepoint) %>% 
    separate(otu, into = c("bactname", "OTUnumber"), sep = "\\ [(]", remove = FALSE) %>% #Add columns to separate bacteria name from OTU number to utilize ggtext so that only bacteria name is italicized
    mutate(otu_name = glue("*{bactname}* ({OTUnumber}")) %>% #Markdown notation so that only bacteria name is italicized
    pull(otu_name)
}
`sig_otu_day-1` <- get_source_otus(-1) 
`sig_otu_day0` <- get_source_otus(0) 
`sig_otu_day1` <- get_source_otus(1)
#OTUs that consistently vary across sources of mice on day-1, 0, or 1
shared_sig_otus_Dn1toD1 <- intersect_all(`sig_otu_day-1`, sig_otu_day0, sig_otu_day1)


#OTUs that were impacted by clindamycin treatment
`sig_otu_pairs` <- read_tsv("data/process/otu_stats_dn1to0.tsv") %>% 
  filter(p.value.adj < 0.05) %>% 
  separate(otu, into = c("bactname", "OTUnumber"), sep = "\\ [(]", remove = FALSE) %>% #Add columns to separate bacteria name from OTU number to utilize ggtext so that only bacteria name is italicized
  mutate(otu_name = glue("*{bactname}* ({OTUnumber}")) %>% #Markdown notation so that only bacteria name is italicized
  pull(otu_name)

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
  
  #Identify otus that overlap with the key_otus that show up across multiple comparisons
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
  source_overlap <- paste(source_df %>% filter(overlap == TRUE) %>% pull(otus), collapse="<br>")
#  source_overlap <- source_df %>% filter(overlap == TRUE) %>% select(otus)
  source_unique <- paste(source_df %>% filter(overlap == FALSE) %>% pull(otus), collapse="<br>")
  clind_overlap <- paste(clind_df %>% filter(overlap == TRUE) %>% pull(otus), collapse="<br>")
  clind_unique <- paste(clind_df %>% filter(overlap == FALSE) %>% pull(otus), collapse="<br>")
#  clind_unique <- clind_df %>% filter(overlap == FALSE) %>% select(otus)
  intersect_overlap <- paste(overlap_df %>% filter(overlap == TRUE) %>% pull(otus), collapse="<br>")
  intersect_unique <- paste(overlap_df %>% filter(overlap == FALSE) %>% pull(otus), collapse="<br>")
  
  otus_venn_plot <- ggplot(df.venn) +
    geom_circle(aes(x0 = x, y0 = y, r = 19, fill = labels), alpha = .1, size = 1, colour = 'grey') +
    coord_fixed() +
    theme_void() +
    theme(legend.position = 'none') +
    scale_fill_manual(values = c('cornflowerblue', 'firebrick',  'gold')) +
    scale_colour_manual(values = c('cornflowerblue', 'firebrick', 'gold'), guide = FALSE) +
    labs(fill = NULL) +
    annotate("text", x = df_venn_otus$x, y = df_venn_otus$y, label = df_venn_otus[,1], size = 5)+
    annotate("text", x = c(-10, 10), y = c(10, 10), label = c("Source", "Clindamycin"), size = 5)+
    annotate(geom='richtext', label = source_overlap, x = -19, y = 0.5, size = 2.8, fill = NA, label.color = NA, color = "firebrick")+
    annotate(geom='richtext', label = source_unique, x = -19, y = -13, size = 2.8, fill = NA, label.color = NA, color = "grey27")+
    annotate(geom='richtext', label = clind_overlap, x = 19, y = 0.5, size = 2.8, fill = NA, label.color = NA, color = "firebrick")+
    annotate(geom='richtext', label = clind_unique, x = 19, y = -13, size = 2.8, fill = NA, label.color = NA, color = "grey27")+
    annotate(geom='richtext', label = intersect_overlap, x = 0, y = -7.5, size = 2.8, fill = NA, label.color = NA, color = "firebrick")+
    annotate(geom='richtext', label = intersect_unique, x = 0, y = -13, size = 2.8, fill = NA, label.color = NA, color = "grey27")+
    annotate("text", x = 0, y = 14, label = title, size = 5)
}

# Venn diagram of Day -1 overlapping OTUs----
dn1_otus_venn_plot <- venn_otus(dayn1_and_interp_dn1_o, paired_and_interp_dn1_o, "Day -1 model OTU comparisons")
save_plot("results/figures/venn_dn1_otus.png", dn1_otus_venn_plot, base_aspect_ratio = 1.8)

# Venn diagram of Day 0 overlapping OTUs----
d0_otus_venn_plot <- venn_otus(day0_and_interp_d0_o, paired_and_interp_d0_o, "Day 0 model OTU comparisons")
save_plot("results/figures/venn_d0_otus.png", d0_otus_venn_plot, base_aspect_ratio = 1.8)

# Venn diagram of Day -1 overlapping OTUs----
d1_otus_venn_plot <- venn_otus(day1_and_interp_d1_o, paired_and_interp_d1_o, "Day 1 model OTU comparisons")
save_plot("results/figures/venn_d1_otus.png", d1_otus_venn_plot, base_aspect_ratio = 1.8)

# Venn diagram of combined overlapping OTUs----
combined_otus_venn_plot <- venn_otus(dayn1to1_and_interp_combined, paired_and_interp_combined, "OTU comparisons for day -1, 0, and 1 models")
save_plot("results/figures/venn_overall_otus.png", combined_otus_venn_plot, base_aspect_ratio = 1.8)

#Comparison of combined Venn diagram OTUs to significant OTUs that varied by source on day -1, 0, and 1:
intersect_all(`shared_sig_otus_Dn1toD1`, `dayn1to1_and_interp_combined`)
intersect_all(`shared_sig_otus_Dn1toD1`, `paired_and_interp_combined`)
