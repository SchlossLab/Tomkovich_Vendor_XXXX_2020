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

#Extract sample ids from distance matrix and join to metadata in order to test relevant variables
all_variables <- tibble(id = attr(all_dist, "Labels")) %>% 
  left_join(metadata, by = "id") 

#Get factor levels for mouse_id variable:
mouse_id_levels <- unique(as.factor(all_data$mouse_id))

variables <- all_data %>% 
  select(-starts_with("Otu")) %>% #get rid of Otu variables
  mutate(experiment=factor(experiment, levels=c("1", "2")), # Make sure variables of interest are treated as factors
         vendor=factor(vendor, levels=c("Schloss", "Young", "Jackson", "Charles River", "Taconic", "Envigo")),
         run=factor(run, levels=c("run_1", "run_2")),
         day=factor(day, levels=c("-1", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9")),
         mouse_id=factor(mouse_id, levels =mouse_id_levels))

#Statistical Analysis----

run_adonis <- adonis(all_dist~ run + mouse_id, data = variables, method = "bray", permutations = 9999)
# P = 1e-04
# P = 1e-04
exp_adonis <- adonis(all_dist~ experiment + mouse_id, data = variables, method = "bray", permutations = 9999)
# P = 1e-04
# P = 1e-04
vendor_adonis <- adonis(all_dist~ vendor + mouse_id, data = variables, method = "bray", permutations = 9999)
# P = 1e-04 
# P = 1e-04
day_adonis <- adonis(all_dist~ day + mouse_id, data = variables, method = "bray", permutations = 9999)
# P = 1e-04
# P = 1e-04

all_adonis <- adonis(all_dist~vendor + mouse_id + day + experiment + run, data = variables, method = "bray", permutations = 9999)

tibble(effects = c("vendor", "mouse_id", "day", "experiment", "run"),
        r_sq = all_adonis$aov.tab$R2[1:5],
        p = all_adonis$aov.tab$Pr[1:5]) %>% 
  write_tsv("data/process/adonis_all.tsv")

#Test ordering (Terms added sequentially: first to last)
all_adonis_test <- adonis(all_dist~experiment + run + vendor + mouse_id + day, data = variables, method = "bray", permutations = 9999)

tibble(effects = c("experiment", "run", "vendor", "mouse_id", "day"),
       r_sq = all_adonis$aov.tab$R2[1:5],
       p = all_adonis$aov.tab$Pr[1:5]) %>% 
  write_tsv("data/process/adonis_all.tsv")
