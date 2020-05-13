source("code/functions.R")

pcoa_data <- read_tsv("data/process/vendors.subsample.thetayc.ave.pcoa.axes") %>%
  select(group, axis1, axis2) %>% #Limit to 2 PCoA axes
  rename(id = group) %>% #group is the same as id in the metadata data frame
  right_join(metadata, by= "id") %>% #merge metadata and PCoA data frames
  filter(!is.na(axis1)) #Remove all samples that weren't sequenced or were sequenced and didn't make the subsampling cutoff


#Determine number of samples with sequence data for each day of the experiment:
seq_data_per_day <- pcoa_data %>% group_by(day) %>% 
  count() %>% arrange(desc(n))
#We only have sequence data for 18 mice on Day 2 and 21 mice on Day 8. 

#Function to plot pcoa data for all vendors----
plot_pcoa <- function(df){
  ggplot(df, aes(x=axis1, y=axis2, color = vendor, alpha = day)) +
    geom_point(size=2) +
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_vendors,
                        labels=color_vendors)+
    scale_alpha_continuous(range = c(.3, 1),
                           breaks= c(2, 4, 6, 8, 10),
                           labels=c(2, 4, 6, 8, 10))+
    coord_fixed() + 
    labs(x="PCoA 1",
         y="PCoA 2",
         color= "Vendor",
         alpha= "Day") +
    theme_classic()
}

#PCoA plot that combines the 2 experiments and save the plot----  
pcoa_plot_combined <- plot_pcoa(pcoa_data) 
save_plot(filename = paste0("results/figures/pcoa.png"), pcoa_plot_combined)

#Statistical Analysis
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

run_adonis <- adonis(all_dist~ run, strata = variables$mouse_id, data = variables, permutations = 9999)
# P = 3e-04
# R2 = .00638
exp_adonis <- adonis(all_dist~ experiment, strata = variables$mouse_id, data = variables, permutations = 9999)
# P = 1 NS
# R2 = .01712
source_adonis <- adonis(all_dist~ source, strata = variables$mouse_id, data = variables, permutations = 9999)
# P = 1 NS
# R2 = .35
day_adonis <- adonis(all_dist~ day, strata = variables$mouse_id, data = variables, permutations = 9999)
# P = 1e-04
# R2 = .11025

all_adonis <- adonis(all_dist~(source/(unique_cage*experiment*run))*day, strata = variables$mouse_id, data = variables, permutations = 9999)

tibble(effects = c("source", "day", "source:unique_cage", "source:run", "source:day", "source:unique_cage:run", "source:unique_cage:day"),
        r_sq = all_adonis$aov.tab$R2[1:7],
        p = all_adonis$aov.tab$Pr[1:7]) %>% 
  write_tsv("data/process/adonis_all.tsv")

#Function to plot pcoa data for all sources of mice at a specific timepoint----
plot_pcoa <- function(df, timepoint){
  plot <- ggplot(df, aes(x=axis1, y=axis2, color = vendor)) +
    geom_point(size=4, alpha = 0.4) +
    scale_colour_manual(name=NULL,
                        values=color_scheme,
                        breaks=color_vendors,
                        labels=color_vendors)+
    coord_fixed() + 
    labs(x="PCoA 1",
         y="PCoA 2",
         color= "Source") +
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5),
          text = element_text(size = 16),
          legend.position = "bottom")
}

# Read in thetayc distance matrix that represents day -1 sequenced samples----
dn1_dist <- read_dist("data/mothur/d-1/vendors.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.subsample.thetayc.0.03.lt.ave.dist")

#Extract sample ids from distance matrix and join to metadata in order to test relevant variables
dn1_variables <- tibble(id = attr(dn1_dist, "Labels")) %>% 
  left_join(variables, by = "id") 

#Unique cage & source are not completely independent so these variables should be nested along with experiment, which cage is nested within
#All of these samples were on the same MiSeq run so that variable does not need to be examined
dn1_adonis <- adonis(dn1_dist~source/(unique_cage*experiment), data = dn1_variables, permutations = 9999)

dn1_results <- tibble(effects = c("source", "source:unique_cage"),
       r_sq = dn1_adonis$aov.tab$R2[1:2],
       p = dn1_adonis$aov.tab$Pr[1:2]) %>% 
  mutate(day = -1)

dn1_pcoa <- read_tsv("data/mothur/d-1/vendors.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.subsample.thetayc.0.03.lt.ave.pcoa.axes") %>%
  select(group, axis1, axis2) %>% #Limit to 2 PCoA axes
  rename(id = group) %>% #group is the same as id in the metadata data frame
  right_join(metadata, by= "id") %>% #merge metadata and PCoA data frames
  filter(!is.na(axis1)) #Remove all samples that weren't sequenced or were sequenced and didn't make the subsampling cutoff

dn1 <- plot_pcoa(dn1_pcoa, -1)+
  ggtitle("Baseline")+ #Title plot
  theme(plot.title = element_text(hjust = 0.5)) #Center plot titile
save_plot(filename = paste0("results/figures/pcoa_day-1.png"), dn1)

# Read in thetayc distance matrix that represents day 0 sequenced samples----
d0_dist <- read_dist("data/mothur/d0/vendors.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.subsample.thetayc.0.03.lt.ave.dist")

#Extract sample ids from distance matrix and join to metadata in order to test relevant variables
d0_variables <- tibble(id = attr(d0_dist, "Labels")) %>% 
  left_join(variables, by = "id") 

#Unique cage & source are not completely independent so these variables should be nested along with experiment, which cage is nested within
#All of these samples were on the same MiSeq run so that variable does not need to be examined
d0_adonis <- adonis(d0_dist~source/(unique_cage*experiment), data = d0_variables, permutations = 9999)

d0_results <- tibble(effects = c("source", "source:unique_cage"),
       r_sq = d0_adonis$aov.tab$R2[1:2],
       p = d0_adonis$aov.tab$Pr[1:2]) %>% 
  mutate(day = 0)

d0_pcoa <- read_tsv("data/mothur/d0/vendors.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.subsample.thetayc.0.03.lt.ave.pcoa.axes") %>%
  select(group, axis1, axis2) %>% #Limit to 2 PCoA axes
  rename(id = group) %>% #group is the same as id in the metadata data frame
  right_join(metadata, by= "id") %>% #merge metadata and PCoA data frames
  filter(!is.na(axis1)) #Remove all samples that weren't sequenced or were sequenced and didn't make the subsampling cutoff

d0 <- plot_pcoa(d0_pcoa, 0)+
  ggtitle("Clindamycin")+ #Title plot
  theme(plot.title = element_text(hjust = 0.5)) #Center plot titile
save_plot(filename = paste0("results/figures/pcoa_day0.png"), d0)


# Read in thetayc distance matrix that represents day 1 sequenced samples----
d1_dist <- read_dist("data/mothur/d1/vendors.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.subsample.thetayc.0.03.lt.ave.dist")

#Extract sample ids from distance matrix and join to metadata in order to test relevant variables
d1_variables <- tibble(id = attr(d1_dist, "Labels")) %>% 
  left_join(variables, by = "id") 

#Unique cage & source are not completely independent so these variables should be nested along with experiment, which cage is nested within
#All of these samples were on the same MiSeq run so that variable does not need to be examined
d1_adonis <- adonis(d1_dist~source/(unique_cage*experiment), data = d1_variables, permutations = 9999)

d1_results <- tibble(effects = c("source", "source:unique_cage"),
       r_sq = d1_adonis$aov.tab$R2[1:2],
       p = d1_adonis$aov.tab$Pr[1:2]) %>% 
  mutate(day = 1)

d1_pcoa <- read_tsv("data/mothur/d1/vendors.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.subsample.thetayc.0.03.lt.ave.pcoa.axes") %>%
  select(group, axis1, axis2) %>% #Limit to 2 PCoA axes
  rename(id = group) %>% #group is the same as id in the metadata data frame
  right_join(metadata, by= "id") %>% #merge metadata and PCoA data frames
  filter(!is.na(axis1)) #Remove all samples that weren't sequenced or were sequenced and didn't make the subsampling cutoff

d1 <- plot_pcoa(d1_pcoa, 1) +
  ggtitle("Post-infection")+ #Title plot
  theme(plot.title = element_text(hjust = 0.5)) #Center plot titile
save_plot(filename = paste0("results/figures/pcoa_day1.png"), d1)


#Make combined table of adonis results for D-1, 0, and 1
rbind(dn1_results, d0_results, d1_results) %>% 
write_tsv("data/process/adonis_dn1-1.tsv")

#Examine initial day -1 vendor communities separately----

#Function to plot pcoa data for each source of mice at day -1, colored according to experiment with shapes to differentiate cage----
pcoa_vendor <- function(df, source, label){
  plot <- ggplot(df, aes(x=axis1, y=axis2, color = experiment, shape = cage)) +
    geom_point(size=4, alpha = 0.4) +
    coord_fixed() + 
    labs(title = label,
         x="PCoA 1",
         y="PCoA 2",
         color= "Experiment",
         shape = "Cage") +
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5),
          text = element_text(size = 16))
  save_plot(filename = paste0("results/figures/pcoa_dn1_", source,".png"), plot)
}

# Read in thetayc distance matrix that represents Day -1 Schloss samples----
s_dn1_dist <- read_dist("data/mothur/d-1/schloss/vendors.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.subsample.thetayc.0.03.lt.ave.dist")

#Extract sample ids from distance matrix and join to metadata in order to test relevant variables
s_dn1_variables <- tibble(id = attr(s_dn1_dist, "Labels")) %>% 
  left_join(variables, by = "id") 

#Unique cage & source are not completely independent so these variables should be nested along with experiment, which cage is nested within
#All of these samples were on the same MiSeq run so that variable does not need to be examined. Generates set of all permutations, which is less than minperm.
s_dn1_adonis <- adonis(s_dn1_dist~experiment*unique_cage, data = s_dn1_variables, complete=TRUE)

s_results <- tibble(effects = c("experiment", "unique_cage"),
       r_sq = s_dn1_adonis$aov.tab$R2[1:2],
       p = s_dn1_adonis$aov.tab$Pr[1:2]) %>% 
  mutate(source= "Schloss")
    
s_dn1_pcoa <- read_tsv("data/mothur/d-1/schloss/vendors.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.subsample.thetayc.0.03.lt.ave.pcoa.axes") %>%
  select(group, axis1, axis2) %>% #Limit to 2 PCoA axes
  rename(id = group) %>% #group is the same as id in the metadata data frame
  right_join(metadata, by= "id") %>% #merge metadata and PCoA data frames
  filter(!is.na(axis1)) #Remove all samples that weren't sequenced or were sequenced and didn't make the subsampling cutoff

s_dn1 <- pcoa_vendor(s_dn1_pcoa, "schloss", "Schloss")

# Read in thetayc distance matrix that represents Day -1 Young samples----
y_dn1_dist <- read_dist("data/mothur/d-1/young/vendors.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.subsample.thetayc.0.03.lt.ave.dist")

#Extract sample ids from distance matrix and join to metadata in order to test relevant variables
y_dn1_variables <- tibble(id = attr(y_dn1_dist, "Labels")) %>% 
  left_join(variables, by = "id") 

#Unique cage & source are not completely independent so these variables should be nested along with experiment, which cage is nested within
#All of these samples were on the same MiSeq run so that variable does not need to be examined. Generates set of all permutations, which is less than minperm.
y_dn1_adonis <- adonis(y_dn1_dist~experiment*unique_cage, data = y_dn1_variables, complete=TRUE)

y_results <- tibble(effects = c("experiment", "unique_cage"),
       r_sq = y_dn1_adonis$aov.tab$R2[1:2],
       p = y_dn1_adonis$aov.tab$Pr[1:2]) %>% 
  mutate(source= "Young")

y_dn1_pcoa <- read_tsv("data/mothur/d-1/young/vendors.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.subsample.thetayc.0.03.lt.ave.pcoa.axes") %>%
  select(group, axis1, axis2) %>% #Limit to 2 PCoA axes
  rename(id = group) %>% #group is the same as id in the metadata data frame
  right_join(metadata, by= "id") %>% #merge metadata and PCoA data frames
  filter(!is.na(axis1)) #Remove all samples that weren't sequenced or were sequenced and didn't make the subsampling cutoff

y_dn1 <- pcoa_vendor(y_dn1_pcoa, "young", "Young")

# Read in thetayc distance matrix that represents Day -1 Jackson samples----
j_dn1_dist <- read_dist("data/mothur/d-1/jackson/vendors.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.subsample.thetayc.0.03.lt.ave.dist")

#Extract sample ids from distance matrix and join to metadata in order to test relevant variables
j_dn1_variables <- tibble(id = attr(j_dn1_dist, "Labels")) %>% 
  left_join(variables, by = "id") 

#Unique cage & source are not completely independent so these variables should be nested along with experiment, which cage is nested within
#All of these samples were on the same MiSeq run so that variable does not need to be examined. Generates set of all permutations, which is less than minperm.
j_dn1_adonis <- adonis(j_dn1_dist~experiment*unique_cage, data = j_dn1_variables, complete=TRUE)

j_results <- tibble(effects = c("experiment", "unique_cage"),
       r_sq = j_dn1_adonis$aov.tab$R2[1:2],
       p = j_dn1_adonis$aov.tab$Pr[1:2]) %>% 
  mutate(source= "Jackson")

j_dn1_pcoa <- read_tsv("data/mothur/d-1/jackson/vendors.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.subsample.thetayc.0.03.lt.ave.pcoa.axes") %>%
  select(group, axis1, axis2) %>% #Limit to 2 PCoA axes
  rename(id = group) %>% #group is the same as id in the metadata data frame
  right_join(metadata, by= "id") %>% #merge metadata and PCoA data frames
  filter(!is.na(axis1)) #Remove all samples that weren't sequenced or were sequenced and didn't make the subsampling cutoff

j_dn1 <- pcoa_vendor(j_dn1_pcoa, "jackson", "Jackson")

# Read in thetayc distance matrix that represents Day -1 Charles River samples----
c_dn1_dist <- read_dist("data/mothur/d-1/charles/vendors.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.subsample.thetayc.0.03.lt.ave.dist")

#Extract sample ids from distance matrix and join to metadata in order to test relevant variables
c_dn1_variables <- tibble(id = attr(c_dn1_dist, "Labels")) %>% 
  left_join(variables, by = "id") 

#Unique cage & source are not completely independent so these variables should be nested along with experiment, which cage is nested within
#All of these samples were on the same MiSeq run so that variable does not need to be examined. Generates set of all permutations, which is less than minperm.
c_dn1_adonis <- adonis(c_dn1_dist~experiment*unique_cage, data = c_dn1_variables, complete=TRUE)

c_results <- tibble(effects = c("experiment", "unique_cage"),
       r_sq = c_dn1_adonis$aov.tab$R2[1:2],
       p = c_dn1_adonis$aov.tab$Pr[1:2]) %>% 
  mutate(source= "Charles River")

c_dn1_pcoa <- read_tsv("data/mothur/d-1/charles/vendors.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.subsample.thetayc.0.03.lt.ave.pcoa.axes") %>%
  select(group, axis1, axis2) %>% #Limit to 2 PCoA axes
  rename(id = group) %>% #group is the same as id in the metadata data frame
  right_join(metadata, by= "id") %>% #merge metadata and PCoA data frames
  filter(!is.na(axis1)) #Remove all samples that weren't sequenced or were sequenced and didn't make the subsampling cutoff

c_dn1 <- pcoa_vendor(c_dn1_pcoa, "charles_river", "Charles River")

# Read in thetayc distance matrix that represents Day -1 Taconic samples----
t_dn1_dist <- read_dist("data/mothur/d-1/taconic/vendors.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.subsample.thetayc.0.03.lt.ave.dist")

#Extract sample ids from distance matrix and join to metadata in order to test relevant variables
t_dn1_variables <- tibble(id = attr(t_dn1_dist, "Labels")) %>% 
  left_join(variables, by = "id") 

#Unique cage & source are not completely independent so these variables should be nested along with experiment, which cage is nested within
#All of these samples were on the same MiSeq run so that variable does not need to be examined. Generates set of all permutations, which is less than minperm.
t_dn1_adonis <- adonis(t_dn1_dist~experiment*unique_cage, data = t_dn1_variables, complete=TRUE)

t_results <- tibble(effects = c("experiment", "unique_cage"),
       r_sq = t_dn1_adonis$aov.tab$R2[1:2],
       p = t_dn1_adonis$aov.tab$Pr[1:2]) %>% 
  mutate(source= "Taconic")

t_dn1_pcoa <- read_tsv("data/mothur/d-1/taconic/vendors.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.subsample.thetayc.0.03.lt.ave.pcoa.axes") %>%
  select(group, axis1, axis2) %>% #Limit to 2 PCoA axes
  rename(id = group) %>% #group is the same as id in the metadata data frame
  right_join(metadata, by= "id") %>% #merge metadata and PCoA data frames
  filter(!is.na(axis1)) #Remove all samples that weren't sequenced or were sequenced and didn't make the subsampling cutoff

t_dn1 <- pcoa_vendor(t_dn1_pcoa, "taconic", "Taconic")

# Read in thetayc distance matrix that represents Day -1 Envigo samples----
e_dn1_dist <- read_dist("data/mothur/d-1/envigo/vendors.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.subsample.thetayc.0.03.lt.ave.dist")

#Extract sample ids from distance matrix and join to metadata in order to test relevant variables
e_dn1_variables <- tibble(id = attr(e_dn1_dist, "Labels")) %>% 
  left_join(variables, by = "id") 

#Unique cage & source are not completely independent so these variables should be nested along with experiment, which cage is nested within
#All of these samples were on the same MiSeq run so that variable does not need to be examined. Generates set of all permutations, which is less than minperm.
e_dn1_adonis <- adonis(e_dn1_dist~experiment*unique_cage, data = e_dn1_variables, complete=TRUE)

e_results <- tibble(effects = c("experiment", "unique_cage"),
       r_sq = e_dn1_adonis$aov.tab$R2[1:2],
       p = e_dn1_adonis$aov.tab$Pr[1:2]) %>% 
  mutate(source= "Envigo")

e_dn1_pcoa <- read_tsv("data/mothur/d-1/envigo/vendors.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.subsample.thetayc.0.03.lt.ave.pcoa.axes") %>%
  select(group, axis1, axis2) %>% #Limit to 2 PCoA axes
  rename(id = group) %>% #group is the same as id in the metadata data frame
  right_join(metadata, by= "id") %>% #merge metadata and PCoA data frames
  filter(!is.na(axis1)) #Remove all samples that weren't sequenced or were sequenced and didn't make the subsampling cutoff

e_dn1 <- pcoa_vendor(e_dn1_pcoa, "envigo", "Envigo")

#Merge adonis results for each source of mice together to create one final results table
rbind(s_results, y_results, j_results, c_results, t_results, e_results) %>% 
write_tsv("data/process/adonis_dn1_source.tsv")

