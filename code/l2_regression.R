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

######################## RUN PIPELINE #############################
# Get the cv and test AUCs for 100 data-splits

 input <- commandArgs(trailingOnly=TRUE) # receive input from model
 # Get variables from command line
 seed <- as.numeric(input[1])
 path <- input[2]
 data <- read.csv(input[3])

 if(!dir.exists(path)){
 	dir.create(path, recursive=TRUE)
 }

 start_time <- Sys.time()

 set.seed(seed)
 get_RMSE_R2_MAE(data, seed, path) # model will be "L2 regularized L2-loss support vector regression"

 end_time <- Sys.time()
 print(end_time - start_time)
###################################################################

# #Test of day -1 input data
# seed <- 0
# path <- "data/process/regression/dayn1"
# data <- read.csv("data/process/regress_input_day-1_data.csv")
# get_RMSE_R2_MAE(data, seed, path) # model will be "L2 regularized L2-loss support vector regression"
# 
# #Test of day 0 input data
# seed <- 0
# path <- "data/process/regression/day0"
# data <- read.csv("data/process/regress_input_day0_data.csv")
# get_RMSE_R2_MAE(data, seed, path) # model will be "L2 regularized L2-loss support vector regression"
# 
# #Test of day 1 input data
# seed <- 0
# path <- "data/process/regression/day1"
# data <- read.csv("data/process/regress_input_day1_data.csv")
# get_RMSE_R2_MAE(data, seed, path) # model will be "L2 regularized L2-loss support vector regression"


