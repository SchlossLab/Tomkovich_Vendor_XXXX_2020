#!/bin/bash

#######################################################################################
# This script will:
#   1. Take single .csv files that have the model result for one datasplit
#   2. Combine them together to have the results for 100 datasplits in one .csv file
#   3. We don't keep the header of each file when combined but only once.

# In the end, the combined_best file must be 101 lines. 1st line is the header and the 100 lines have the data of 100 files.
#             the combined_all file must have 100*(hyper-parameter number)+1 lines.
########################################################################################


SEARCH_DIR=data/process/regression/day1
FINAL_DIR=data/process/regression
input=day1 #Input Otu data that was used in model to predict C. difficile colonization status on day7

# Keep the first line of File1 and remove the first line of all the others and combine


head -1 $SEARCH_DIR/all_cost.1.csv  > $SEARCH_DIR/combined_all_cost_"$input".csv; tail -n +2 -q $SEARCH_DIR/all_cost.*.csv >> $SEARCH_DIR/combined_all_cost_"$input".csv
head -1 $SEARCH_DIR/optimum_cost.1.csv  > $SEARCH_DIR/combined_optimum_cost_"$input".csv; tail -n +2 -q $SEARCH_DIR//optimum_cost.*.csv >> $SEARCH_DIR/combined_optimum_cost_"$input".csv


mv $SEARCH_DIR/combined_all_cost_"$input".csv $FINAL_DIR/combined_all_cost_"$input".csv
mv $SEARCH_DIR/combined_optimum_cost_"$input".csv $FINAL_DIR/combined_optimum_cost__"$input".csv
