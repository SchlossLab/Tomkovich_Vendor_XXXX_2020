# Regression to predict C. diff cfu on D5 using Otu data from day 1
#!/bin/bash
# SEED for full analysis
#for S in {0..99}
for S in {0..99}
do
    Rscript code/l2_regression.R $S "data/process/regression/day1" "data/process/regress_input_day1_data.csv"
done
