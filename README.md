Download the [latest release](https://github.com/SchlossLab/new_project/releases/latest) to the directory and decompress


## Initial gut microbiota influences Clostridioides difficile clearance after clindamycin perturbation in mice

ABSTRAAAAAAAAAAAAACT




### Overview

	project
	|- README          # the top level description of content (this doc)
	|- CONTRIBUTING    # instructions for how to contribute to your project
	|- LICENSE         # the license for this project
	|
	|- submission/
	| |- study.Rmd    # executable Rmarkdown for this study, if applicable
	| |- study.md     # Markdown (GitHub) version of the *.Rmd file
	| |- study.tex    # TeX version of *.Rmd file
	| |- study.pdf    # PDF version of *.Rmd file
	| |- header.tex   # LaTeX header file to format pdf version of manuscript
	| |- references.bib # BibTeX formatted references
	| |- XXXX.csl     # csl file to format references for journal XXX
	|
	|- data           # raw and primary data, are not changed once created
	| |- references/  # reference files to be used in analysis
	| |- raw/         # raw data, will not be altered
	| |- mothur/      # mothur processed data
	| +- process/     # cleaned data, will not be altered once created;
	|                 # will be committed to repo
	|
	|- code/          # any programmatic code
	|
	|- results        # all output from workflows and analyses
	| |- tables/      # text version of tables to be rendered with kable in R
	| |- figures/     # graphs, likely designated for manuscript figures
	| +- pictures/    # diagrams, images, and other non-graph graphics
	|
	|- exploratory/   # exploratory data analysis for study
	| |- notebook/    # preliminary analyses
	| +- scratch/     # temporary files that can be safely deleted or lost
	|
	+- Makefile       # executable Makefile for this study, if applicable


### How to regenerate this repository

#### Dependencies and locations
* Gnu Make should be located in the user's PATH
* mothur (v1.43.0) should be located in the user's PATH
* R (v. 3.5.2) should be located in the user's PATH
* R packages:
    * broom v0.5.0
    * tidyverse_1.3.0
    * cowplot v0.9.3
    * magick v2.0
    * vegan v2.5-2
    * reshape2 v1.4.3
    * knitr v1.20
    * rmarkdown v1.10
    * gtools v3.8.1
    * ggpubr v.0.2.4
    *ggforce v0.1.3
* Analysis assumes the use of 8 processors  


#### Running analysis
Download 16S rRNA sequencing dataset from the NCBI Sequence Read Archive (BioProject Accession no. PRJNA608529).
```
git clone https://github.com/SchlossLab/Tomkovich_vendor_difs_XXXX_2020
make write.paper
```
Transfer 16S rRNA sequencing fastq.gz files into Tomkovich_vendor_difs_XXX_2020/data/raw
```
cd Tomkovich_vendor_difs_XXXX_2020
```
Obtain the SILVA reference alignment from version 132 described at http://blog.mothur.org/2018/01/10/SILVA-v132-reference-files/. We will use the SEED v. 132, which contain 12,083 bacterial sequences. This also contains the reference taxonomy. We will limit the databases to only include bacterial sequences.
```
wget -N https://mothur.org/w/images/7/71/Silva.seed_v132.tgz
tar xvzf Silva.seed_v132.tgz silva.seed_v132.align silva.seed_v132.tax
mothur "#get.lineage(fasta=silva.seed_v132.align, taxonomy=silva.seed_v132.tax, taxon=Bacteria);degap.seqs(fasta=silva.seed_v132.pick.align, processors=8)"
mv silva.seed_v132.pick.align data/references/silva.seed.align
rm Silva.seed_v132.tgz silva.seed_v132.*
#Narrow to v4 region
mothur "#pcr.seqs(fasta=data/references/silva.seed.align, start=11894, end=25319, keepdots=F, processors=8)"
mv data/references/silva.seed.pcr.align data/references/silva.v4.align
```
Obtain the RDP reference taxonomy. The current version is v11.5 and we use a "special" pds version of the database files, which are described at http://blog.mothur.org/2017/03/15/RDP-v16-reference_files/.
```
wget -N https://www.mothur.org/w/images/c/c3/Trainset16_022016.pds.tgz
tar xvzf Trainset16_022016.pds.tgz trainset16_022016.pds
mv trainset16_022016.pds/* data/references/
rm -rf trainset16_022016.pds
rm Trainset16_022016.pds.tgz
```
Obtain the Zymo mock community data; note that Zymo named the 5 operon of Salmonella twice instead of the 7 operon.
```
wget -N https://s3.amazonaws.com/zymo-files/BioPool/ZymoBIOMICS.STD.refseq.v2.zip
unzip ZymoBIOMICS.STD.refseq.v2.zip
rm ZymoBIOMICS.STD.refseq.v2/ssrRNAs/*itochondria_ssrRNA.fasta #V4 primers don't come close to annealing to these
cat ZymoBIOMICS.STD.refseq.v2/ssrRNAs/*fasta > zymo_temp.fasta
sed '0,/Salmonella_enterica_16S_5/{s/Salmonella_enterica_16S_5/Salmonella_enterica_16S_7/}' zymo_temp.fasta > zymo.fasta
mothur "#align.seqs(fasta=zymo.fasta, reference=data/references/silva.v4.align, processors=12)"
mv zymo.align data/references/zymo_mock.align
rm -rf zymo* ZymoBIOMICS.STD.refseq.v2* zymo_temp.fasta
```
Run the 16S rRNA sequnce data through mothur.
Create vendors.files:
```
mothur code/make_file_fix_ids.batch
```
Fix the sample ids in vendor.files before proceeding with the rest of the analysis by doing the following:
Make a copy of data/raw/vendors.files and save as a .csv file (data/raw/vendors.files.csv)
Then run the following rscript to fix the ids:
```
Rscript code/fix_ids.R
```
Open up data/raw/vendors.files.fixed.csv and data/raw/vendors.files. Paste rhe fix_id column (select all 439 samples, starting at row 2) from data/raw/vendors.files.fixed.csv over the 1st column with dashes in the group name in data/raw/vendors.files and save data/raw/vendors.files. Now that the sample ids have been corrected, the analysis can proceed.
```
mothur code/get_good_seqs.batch
mothur code/get_error.batch
mothur code/get_shared_otus.batch
mothur code/get.otu.batch
```
Copy files generated by mothur analysis and place into data/process while shortening file names at the same time.
```
bash rename_mothur_outputs
```

Used code/split_analysis_by_day.R to generate lists of ids specific to each sample days and pasted them into the following scripts to create a bray-curtis distance matrix and PCoA specific for a single timepoint.
```
mothur code/d-1_dist_PCoA
mothur code/d0_dist_PCoA
```

For L2 Logistic regression analysis, use the following script to generate the input data:
```

Rscript code/l2_classification_input_data.R
```
Use ML_pipeline_microbiome to perform L2 Logistic regression analysis. Modify model_pipeline.R to specify outcomes. Paste input_*data.csv files into ML_pipeline_microbiome/test/data. Run the pipeline, merge files in between, otherwise they will be overwritten in temp folder:
```
mkdir data/process/day0 data/process/day1 data/process/dayn1
bash code/bash/L2_log_Regression_d0.sh
bash code/bash/d0_cat_csv_files.sh
bash code/bash/L2_log_Regression_d1.sh
code/bash/d1_cat_csv_files.sh
bash code/bash/L2_log_Regression_dn1.sh
code/bash/dn1_cat_csv_files.sh
```
Copy outputs from Logistic Regression classification analysis (5 files in each of the following subfolders data/process/day0, data/process/day1, data/process/dayn1) and paste into data/process/classification.

Rename classification output files to indicate which day of the experiment Otu relative abundances were used in the model to predict C. difficile colonization status on day 7.
```
bash code/rename_classification_outputs
```
Make plots comparing cross-validation and test AUROCs for Logistic Regression models generated from 3 different types of input data (OTUs from Day -1, 0, and 1) and make plots to interpret which OTUs are contributing to the 3 models.
```
Rscript code/class._regress._analysis.R
Rscript code/class_interpretation.R
```

```
make write.paper
```
