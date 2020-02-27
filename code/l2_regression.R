# Load needed libraries
library("tidyverse")
library("caret")
library("pROC")
library("LiblineaR")
library("kernlab")
library("doParallel")
library("gtools")
library("data.table")

# Function to run Begum's pipeline
pipeline <- function(dataset){

# Create vectors to save cv and test AUC values for every data-split
results_total <-  data.frame()
test_aucs <- c()
cv_aucs <- c()

# We are doing the pre-processing to the full dataset and then splitting 80-20
# Scale all features between 0-1
#preProcValues <- preProcess(dataset, method = "range")
#dataTransformed <- predict(preProcValues, dataset)
#Will confirm, but I think Begum said the Preprocessing step wasn't necessary

# remove columns that only appear within one or fewer samples of the training set. These are
# likely to be all zero and will not enter into the model
frequent <- names(which(apply(dataset[, -1] > 0, 2, sum) > 1))

dataTransformed <- dataset %>% select(regress, frequent)

# Do the 80-20 data-split
# Stratified data partitioning %80 training - %20 testing

inTraining <- createDataPartition(dataTransformed$regress, p = .80, list = FALSE)
trainTransformed <- dataTransformed[ inTraining,]
testTransformed  <- dataTransformed[-inTraining,]

#	n_features <- ncol(trainTransformed) - 1
#	if(n_features > 20000) n_features <- 20000

#	if(n_features < 19){ mtry <- 1:6
#	} else { mtry <- floor(seq(1, n_features/3, length=6)) }

#	mtry <- mtry[mtry <= n_features]
  

# cv index to make sure the internal 5-folds are stratified for diagnosis classes and also resampled 100 times.
# 100 repeat internally is necessary to get robust readings of hyperparameter setting performance
  folds <- 5
	cvIndex <- createMultiFolds(factor(trainTransformed$regress), folds, times=100) #returnTrain = T default for multifolds
  cv <- trainControl(method="repeatedcv",
                     number= folds,
                     index = cvIndex,
                     returnResamp="final",
                     classProbs=FALSE,
                     indexFinal=NULL,
                     savePredictions = TRUE)

  # Hyper-parameter tuning budget 
  hyperparameters <- read.csv('data/default_hyperparameters.csv', stringsAsFactors = F) #Source: From https://github.com/SchlossLab/ML_pipeline_microbiome/tree/master/data and narrowed to L2_Linear_SVM
  hyperparameters <- split(hyperparameters$val, hyperparameters$param)
	grid <-  expand.grid(cost = hyperparameters$cost,
	                     Loss = "L2")	# This chooses type=0 for liblinear R package *Need to adapt for L2-regularized L2-loss support vector regression (primal). Type 11 in LiblineaR package
	                     # https://www.rdocumentation.org/packages/LiblineaR/versions/2.10-8/topics/LiblineaR
  method <- "svmLinear3"
  
  # Train model for 1 data-split but with 10fold-10repeat CV
  trained_model <-  train(regress ~ .,
                          data=trainTransformed,
                          method = "svmLinear3",
                          trControl = cv,
                          metric = "RMSE",
                          tuneGrid = grid,
                          family = "binomial",
                          svr_eps = 0.01 #default is 0.1 if not specified https://www.rdocumentation.org/packages/LiblineaR/versions/2.10-8/topics/LiblineaR
  )

  # RMSE value over repeats of the best cost parameter during training
	cv_best <- getTrainPerf(trained_model)

  # Training results for all the cost parameters
  cv_results <- trained_model$results

  # Predictions with the best cost parameter for 1 data-split
  predictions <- predict(trained_model, testTransformed, type = "raw")

  # quality values for test data from 1 datasplit
	test_results <- postResample(testTransformed$regress, predictions)

  results <- list(cv_best=cv_best, cv_all=cv_results, test=test_results)
  return(results)
}

# Function to save the RMSE values and save them as .csv
get_RMSE_R2_MAE <- function(dataset, split_number, dir){

  # Save results of the modeling pipeline as a list
  results <- pipeline(dataset)

  # ------------------------------------------------------------------
  # Create a matrix with cv_aucs and test_aucs from 100 data splits
	colnames(results$cv_best) <- str_replace(colnames(results$cv_best), "Train", "train_")
	names(results$test) <- paste0("test_", names(results$test))


  # Convert to dataframe and add a column noting the model name
  enframe(unlist(c(results$cv_best[1, -4], results$test))) %>%
		spread(name, value) %>%
    write_csv(path=paste0(dir, "/optimum_cost.", split_number, ".csv"))

  # ------------------------------------------------------------------

  # ------------------------------------------------------------------
  # Save all tunes from 100 data splits and corresponding AUCs
  results$cv_all %>%
    write_csv(path=paste0(dir, "/all_cost.", split_number, ".csv"))
  # ------------------------------------------------------------------

}

##############################################################################################
############### Run the actual programs to get the data ######################################
##############################################################################################

#Code for formatting data for running pipeline for the SCFA paper
# read_shared <- function(shared_file_name, min_samples=2, min_abundance=1){
# 
# 	data <- fread(shared_file_name, header=T, colClasses=c(Group="character"))[, -c(1,3)]
# 
# 	if(min_abundance != 1){
# 		abundant <- names(which(apply(data[,-1], 2, sum) >= min_abundance))
# 		data <- select(data, Group, abundant)
# 	}
# 
# 	frequent <- names(which(apply(data[, -1] > 0, 2, sum) >= min_samples))
# 
# 	select(data, Group, frequent)
# 
# }
# 
# get_data <- function(path){
# 
# 	parse_path <- unlist(str_split(unlist(str_split(path, '/'))[3], "_"))
# 	target_scfa <- parse_path[1]
# 	feature_sources <- parse_path[-1]
# 
# 	tax_levels <- c("kingdom", "phylum", "class", "order", "family", "genus")
# 
# 	# Read in metadata
# 	if(any(feature_sources %in% metagenomics)){
# 		data <- read_tsv('data/metadata/zackular_metadata.tsv') %>%
# 						mutate(sample=str_replace(sample, "(\\d{1,2})", " \\1")) %>%
# 						separate(sample, into=c("disease", "subject")) %>%
# 						mutate(disease=str_replace(disease, "^(.).*", "\\1"),
# 										dx=tolower(dx)) %>%
# 						unite(sample, disease, subject, sep="") %>%
# 						select(sample, fit_result, dx)
# 	} else {
# 		data <- read_csv('data/metadata/cross_section.csv', col_types=cols(sample=col_character()))
# 	}
# 
# 	# Read in OTU table and remove label and numOtus columns
# 	if("otu" %in% feature_sources){
# 
# 		data <- read_shared('data/mothur/crc.otu.shared') %>%
# 			inner_join(data, ., by=c("sample"="Group"))
# 
# 	}
# 
# 	if(any(feature_sources %in% tax_levels)){
# 		taxon <- feature_sources[which(feature_sources %in% tax_levels)]
# 		shared_taxon_file <- paste0('data/phylotype/crc.', taxon, ".shared")
# 
# 		data <- read_shared(shared_taxon_file) %>%
# 			inner_join(data, ., by=c("sample"="Group"))
# 
# 	}
# 
# 	# Read in SCFAs spread columns
# 	read_tsv('data/scfa/scfa_composite.tsv', col_types=cols(study_id=col_character())) %>%
# 		spread(key=scfa, value=mmol_kg) %>%
# 		select(study_id, target_scfa) %>%
# 		inner_join(., data, by=c("study_id" = "sample")) %>%
# 		rename(regress = target_scfa) %>%
# 		drop_na() %>%
# 		select(-study_id)
# }

#Begum's suggestion for formatting input data for ML_pipeline_microbiome----
######################## DATA PREPARATION #############################
# Features: 16S rRNA gene sequences(OTUs) in the stool
# Labels: - C. difficile CFU levels in the stool 5 days post-challenge

# Read in metadata and select only sample id, day, and cfu_d5 columns
source("code/functions.R")
metadata <- metadata %>% select(id, day, cfu_d5)
# Read in OTU table and remove label and numOtus columns
shared <- read.delim('data/process/vendors.subsample.shared', header=T, sep='\t') %>%
  select(-label, -numOtus)
# Merge metadata and OTU table.
# Filter merged dataframe to just timepoints we want to use community structure to predict C. difficile CFU 5 days post-challenge
# Then remove the sample ID and day columns
select_timepoint <- function(timepoint){
  data <- inner_join(metadata, shared, by=c("id"="Group")) %>% 
    filter(day == timepoint) %>% 
    select(-id, -day) %>%
    drop_na() %>%
    select(cfu_d5, everything()) %>%
    rename(regress=cfu_d5) %>% #rename column, the cfu_d5 values are what we want to predict
    write_csv(paste0("data/process/regress_input_day", timepoint, "_data.csv"))
}
data_dn1 <- select_timepoint(-1)
data_d0 <- select_timepoint(0)
data_d1 <- select_timepoint(1)

######################## RUN PIPELINE #############################
# Get the cv and test AUCs for 100 data-splits

# if running from make:
# input <- commandArgs(trailingOnly=TRUE) # receive input from model
# # Get variables from command line
# seed <- as.numeric(input[1])
# path <- input[2]
# 
# if(!dir.exists(path)){
# 	dir.create(path, recursive=TRUE)
# }
# 
# start_time <- Sys.time()
# 
# set.seed(seed)
# data <- read.csv(path)
# get_RMSE_R2_MAE(data, seed, path) # model will be "L2 regularized L2-loss support vector regression"
# 
# end_time <- Sys.time()
# print(end_time - start_time)
###################################################################

#Test of day -1 input data
seed <- 0
path <- "data/process/regression/dayn1"
data <- read.csv("data/process/regress_input_day-1_data.csv")
get_RMSE_R2_MAE(data, seed, path) # model will be "L2 regularized L2-loss support vector regression"

#Test of day 0 input data
seed <- 0
path <- "data/process/regression/day0"
data <- read.csv("data/process/regress_input_day0_data.csv")
get_RMSE_R2_MAE(data, seed, path) # model will be "L2 regularized L2-loss support vector regression"

#Test of day 1 input data
seed <- 0
path <- "data/process/regression/day1"
data <- read.csv("data/process/regress_input_day1_data.csv")
get_RMSE_R2_MAE(data, seed, path) # model will be "L2 regularized L2-loss support vector regression"


