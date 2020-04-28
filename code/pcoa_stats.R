source("code/functions.R")

set.seed(19881117) #Match seed used in mothur analysis scripts

#Function to format distance matrix generated with mothur for use in R. 
#Source: Sze et al. mSphere 2019 https://github.com/SchlossLab/Sze_PCRSeqEffects_mSphere_2019/blob/master/code/vegan_analysis.R
read_dist <- function(dist_file_name){
  
  linear_data <- scan(dist_file_name, what="character", sep="\n", quiet=TRUE)
  
  n_samples <- as.numeric(linear_data[1])
  linear_data <- linear_data[-1]
  
  samples <- str_replace(linear_data, "\t.*", "")
  linear_data <- str_replace(linear_data, "[^\t]*\t", "")
  linear_data <- linear_data[-1]
  
  distance_matrix <- matrix(0, nrow=n_samples, ncol=n_samples)
  
  for(i in 1:(n_samples-1)){
    row <- as.numeric(unlist(str_split(linear_data[i], "\t")))
    distance_matrix[i+1,1:length(row)] <- row
  }
  
  distance_matrix <- distance_matrix + t(distance_matrix)
  rownames(distance_matrix) <- samples
  
  as.dist(distance_matrix)
}

# Read in thetayc distance matrix that represents all sequenced samples----
all_dist <- read_dist("data/process/vendors.subsample.thetayc.ave.dist")

#Get factor levels for mouse_id variable:
mouse_id_levels <- unique(as.factor(metadata$mouse_id))
#49 levels

#Get factor levels for unique_cage variable:
unique_cage_levels <- unique(as.factor(metadata$unique_cage))
#24 levels

#Extract sample ids from distance matrix and join to metadata in order to test relevant variables
all_variables <- tibble(id = attr(all_dist, "Labels")) %>% 
  left_join(metadata, by = "id") %>% 
  # Make sure variables of interest are treated as factors
  mutate(experiment=factor(experiment, levels=c("1", "2")), 
         vendor=factor(vendor, levels=c("Schloss", "Young", "Jackson", "Charles River", "Taconic", "Envigo")),
         run=factor(run, levels=c("run_1", "run_2")),
         day=factor(day, levels=c("-1", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9")),
         mouse_id=factor(mouse_id, levels =mouse_id_levels),  
         unique_cage=factor(unique_cage, levels=unique_cage_levels)) %>% 
         rename(source = vendor) #Rename vendor variable to reflect language used in paper: mice from different sources
  

variables <- all_variables %>% 
  select(-starts_with("Otu")) #get rid of Otu variables
  
#Statistical Analysis----

run_adonis <- adonis(all_dist~ run, strata = variables$mouse_id, data = variables, method = "bray", permutations = 9999)
# P = 3e-04
# R2 = .00638
exp_adonis <- adonis(all_dist~ experiment, strata = variables$mouse_id, data = variables, method = "bray", permutations = 9999)
# P = 1 NS
# R2 = .01712
source_adonis <- adonis(all_dist~ source, strata = variables$mouse_id, data = variables, method = "bray", permutations = 9999)
# P = 1 NS
# R2 = .35
day_adonis <- adonis(all_dist~ day, strata = variables$mouse_id, data = variables, method = "bray", permutations = 9999)
# P = 1e-04
# R2 = .11025

all_adonis <- adonis(all_dist~(source/(unique_cage*experiment*run))*day, strata = variables$mouse_id, data = variables, method = "bray", permutations = 9999)

tibble(effects = c("source", "day", "source:unique_cage", "source:run", "source:day", "source:unique_cage:run", "source:unique_cage:day"),
        r_sq = all_adonis$aov.tab$R2[1:7],
        p = all_adonis$aov.tab$Pr[1:7]) %>% 
  write_tsv("data/process/adonis_all.tsv")

#Function to plot pcoa data for all sources of mice at a specific timepoint----
plot_pcoa <- function(df, timepoint){
  plot <- ggplot(df, aes(x=axis1, y=axis2, color = vendor)) +
    geom_point(size=2, alpha = 0.4) +
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_vendors,
                        labels=color_vendors)+
    coord_fixed() + 
    labs(x="PCoA 1",
         y="PCoA 2",
         color= "Source",
         alpha= "Day") +
    theme_classic()
  save_plot(filename = paste0("results/figures/pcoa_day", timepoint,".png"), plot)
}

# Read in thetayc distance matrix that represents day -1 sequenced samples----
dn1_dist <- read_dist("data/mothur/d-1/vendors.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.subsample.thetayc.0.03.lt.ave.dist")

#Extract sample ids from distance matrix and join to metadata in order to test relevant variables
dn1_variables <- tibble(id = attr(dn1_dist, "Labels")) %>% 
  left_join(variables, by = "id") 

#Unique cage & source are not completely independent so these variables should be nested along with experiment, which cage is nested within
#All of these samples were on the same MiSeq run so that variable does not need to be examined
dn1_adonis <- adonis(dn1_dist~source/(unique_cage*experiment), data = dn1_variables, method = "bray", permutations = 9999)

tibble(effects = c("source", "source:unique_cage"),
       r_sq = dn1_adonis$aov.tab$R2[1:2],
       p = dn1_adonis$aov.tab$Pr[1:2]) %>% 
  write_tsv("data/process/adonis_dn1.tsv")

dn1_pcoa <- read_tsv("data/mothur/d-1/vendors.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.subsample.thetayc.0.03.lt.ave.pcoa.axes") %>%
  select(group, axis1, axis2) %>% #Limit to 2 PCoA axes
  rename(id = group) %>% #group is the same as id in the metadata data frame
  right_join(metadata, by= "id") %>% #merge metadata and PCoA data frames
  filter(!is.na(axis1)) #Remove all samples that weren't sequenced or were sequenced and didn't make the subsampling cutoff

dn1 <- plot_pcoa(dn1_pcoa, -1)

# Read in thetayc distance matrix that represents day 0 sequenced samples----
d0_dist <- read_dist("data/mothur/d0/vendors.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.subsample.thetayc.0.03.lt.ave.dist")

#Extract sample ids from distance matrix and join to metadata in order to test relevant variables
d0_variables <- tibble(id = attr(d0_dist, "Labels")) %>% 
  left_join(variables, by = "id") 

#Unique cage & source are not completely independent so these variables should be nested along with experiment, which cage is nested within
#All of these samples were on the same MiSeq run so that variable does not need to be examined
d0_adonis <- adonis(d0_dist~source/(unique_cage*experiment), data = d0_variables, method = "bray", permutations = 9999)

tibble(effects = c("source", "source:unique_cage"),
       r_sq = d0_adonis$aov.tab$R2[1:2],
       p = d0_adonis$aov.tab$Pr[1:2]) %>% 
  write_tsv("data/process/adonis_d0.tsv")

d0_pcoa <- read_tsv("data/mothur/d0/vendors.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.subsample.thetayc.0.03.lt.ave.pcoa.axes") %>%
  select(group, axis1, axis2) %>% #Limit to 2 PCoA axes
  rename(id = group) %>% #group is the same as id in the metadata data frame
  right_join(metadata, by= "id") %>% #merge metadata and PCoA data frames
  filter(!is.na(axis1)) #Remove all samples that weren't sequenced or were sequenced and didn't make the subsampling cutoff

d0 <- plot_pcoa(d0_pcoa, 0)

# Read in thetayc distance matrix that represents day 1 sequenced samples----
d1_dist <- read_dist("data/mothur/d1/vendors.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.subsample.thetayc.0.03.lt.ave.dist")

#Extract sample ids from distance matrix and join to metadata in order to test relevant variables
d1_variables <- tibble(id = attr(d1_dist, "Labels")) %>% 
  left_join(variables, by = "id") 

#Unique cage & source are not completely independent so these variables should be nested along with experiment, which cage is nested within
#All of these samples were on the same MiSeq run so that variable does not need to be examined
d1_adonis <- adonis(d1_dist~source/(unique_cage*experiment), data = d1_variables, method = "bray", permutations = 9999)

tibble(effects = c("source", "source:unique_cage"),
       r_sq = d1_adonis$aov.tab$R2[1:2],
       p = d1_adonis$aov.tab$Pr[1:2]) %>% 
  write_tsv("data/process/adonis_d1.tsv")

d1_pcoa <- read_tsv("data/mothur/d1/vendors.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.subsample.thetayc.0.03.lt.ave.pcoa.axes") %>%
  select(group, axis1, axis2) %>% #Limit to 2 PCoA axes
  rename(id = group) %>% #group is the same as id in the metadata data frame
  right_join(metadata, by= "id") %>% #merge metadata and PCoA data frames
  filter(!is.na(axis1)) #Remove all samples that weren't sequenced or were sequenced and didn't make the subsampling cutoff

d1 <- plot_pcoa(d1_pcoa, 1)
