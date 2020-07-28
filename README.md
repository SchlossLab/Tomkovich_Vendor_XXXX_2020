Download the [latest release](https://github.com/SchlossLab/new_project/releases/latest) to the directory and decompress


## The initial gut microbiota and response to antibiotic perturbation influence *Clostridioides difficile* colonization in mice

The gut microbiota has a key role in determining susceptibility to *Clostridioides difficile* infections (CDIs). However, much of the mechanistic work examining CDIs in mouse models use animals obtained from a single university colony or vendor. We treated mice from 6 different sources (2 University of Michigan colonies and 4 vendors) with a single clindamycin dose, followed by a *C. difficile* challenge 1 day later and then measured *C. difficile* colonization levels through 9 days post-infection. The microbiota were profiled via 16S rRNA gene sequencing to examine the variation across sources and alterations due to clindamycin treatment and *C. difficile* challenge. While all sources of mice were colonized 1-day post-infection, variation emerged from days 3-7 post-infection with animals from some sources colonized with *C. difficile* for longer and at higher levels. We identified bacteria that varied in relative abundance across sources and throughout the experiment. Some bacteria were consistently impacted by clindamycin treatment in all sources of mice including *Lachnospiraceae*, *Ruminococcaceae*, and *Enterobacteriaceae*. To identify bacteria that were most important to colonization regardless of the source, we created logistic regression models that successfully classified mice based on whether they cleared *C. difficile* by 7 days post-infection using baseline, post-clindamycin, and post-infection community composition data. With these models, we identified 4 bacteria that varied across sources (*Bacteroides*), were altered by clindamycin (*Porphyromonadaceae*), or both (*Enterobacteriaceae* and *Enterococcus*). Microbiota variation across sources better emulates human interindividual variation and can help identify bacterial drivers of phenotypic variation in the context of CDIs.


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
* FFmpeg should be located in the user's PATH
* R (v. 4.0.2) should be located in the user's PATH
* R packages:
    * broom v0.7.0
    * tidyverse_1.3.0
    * cowplot v1.0.0
    * magick v2.4.0
    * vegan v2.5-6
    * reshape2 v1.4.4
    * knitr v1.29
    * rmarkdown v2.3
    * gtools v3.8.2
    * ggpubr v.0.4.0
    * ggforce v0.3.2
    * gganimate v1.0.6
    * writexl v1.3
    * glue v1.4.1
    * ggtext v0.1.0
* Analysis assumes the use of 8 processors  


#### Running analysis
Download 16S rRNA sequencing dataset from the NCBI Sequence Read Archive (BioProject Accession no. PRJNA608529).
```
git clone https://github.com/SchlossLab/Tomkovich_Vendor_XXXX_2020
```
Transfer 16S rRNA sequencing fastq.gz files into Tomkovich_Vendor_XXXX_2020/data/raw
```
cd Tomkovich_Vendor_XXXX_2020
```
Obtain the SILVA reference alignment from version 132 described at https://mothur.org/blog/2018/SILVA-v132-reference-files/. We will use the SEED v. 132, which contain 12,083 bacterial sequences. This also contains the reference taxonomy. We will limit the databases to only include bacterial sequences.
```
wget -N https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.seed_v132.tgz
tar xvzf Silva.seed_v132.tgz silva.seed_v132.align silva.seed_v132.tax
mothur "#get.lineage(fasta=silva.seed_v132.align, taxonomy=silva.seed_v132.tax, taxon=Bacteria);degap.seqs(fasta=silva.seed_v132.pick.align, processors=8)"
mv silva.seed_v132.pick.align data/references/silva.seed.align
rm Silva.seed_v132.tgz silva.seed_v132.*
#Narrow to v4 region
mothur "#pcr.seqs(fasta=data/references/silva.seed.align, start=11894, end=25319, keepdots=F, processors=8)"
mv data/references/silva.seed.pcr.align data/references/silva.v4.align
```
Obtain the RDP reference taxonomy. The current version is v11.5 and we use a "special" pds version of the database files, which are described at https://mothur.org/blog/2017/RDP-v16-reference_files/.
```
wget -N https://mothur.s3.us-east-2.amazonaws.com/wiki/trainset16_022016.pds.tgz
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
Perform *C. difficile* CFU and mouse weight over time analysis and create plots.
```
Rscript code/cfu_plot.R
Rscript code/weight_plot.R
```
Perform alpha diversity analysis and create plots.
```
Rscript code/diversity.R
```
Use code/subset_analysis.R to generate lists of IDs specific to each sample day and pasted them into the following scripts to create a bray-curtis distance matrix and PCoA specific for a single timepoint. Also used to create PCoAs of baseline communities of each mouse colony source.
```
mothur code/d-1_dist_PCoA
mothur code/d0_dist_PCoA
mothur code/d1_dist_PCoA
mothur code/d-1_vendors_dist_PCoA
```
Perform PCoA analysis and create plots. Compare relative distances within and between vendors and plot data.
```
Rscript code/pcoa_stats.R
Rscript code/relative_distance_comparisons.R
```
Perform taxonomic analysis and create plots.
```
Rscript code/taxa_stats.R
```
For L2 Logistic regression analysis, use the following script to generate the input data at the OTU level:
```
Rscript code/l2_classification_input_data.R
```
Use [ML_pipeline_microbiome repository](https://github.com/SchlossLab/ML_pipeline_microbiome) to perform L2 Logistic regression analysis. Modify model_pipeline.R to specify outcomes and a 60:40 data split for cross-validation and testing steps. Acccess [modified version of repository here](https://github.com/tomkoset/ML_pipeline_microbiome). Move Tomkovich_Vendor_XXXX_2020/data/process/classification_input_*data.csv files into ML_pipeline_microbiome/test/data. Run the pipeline as arrayed jobs using batch scripts formatted for Slurm (or whatever scheduler your HPC uses), merge files in between with bash scripts, otherwise they will be overwritten in temp folder:
```
mkdir data/process/dayn1 data/process/day0 data/process/day1
sbatch code/slurm/L2_log_Regression_dn1.sh
bash code/bash/dn1_cat_csv_files.sh
sbatch code/slurm/L2_log_Regression_d0.sh
bash code/bash/d0_cat_csv_files.sh
sbatch code/slurm/L2_log_Regression_d1.sh
bash code/bash/d1_cat_csv_files.sh

```
Copy the 6 folders containing outputs from Logistic Regression classification analysis (there should be 5 files in each folder) and paste into Tomkovich_Vendor_XXXX_2020/data/process/classification/
ML_pipeline_microbiome/data/process/dayn1
ML_pipeline_microbiome/data/process/day0
ML_pipeline_microbiome/data/process/day1

Rename classification output files to indicate which day of the experiment relative abundances were used to train the models to predict *C. difficile* colonization status on day 7.
```
bash code/rename_classification_outputs_60
```
Make plots comparing cross-validation and test AUROCs for classification models generated from 3 different types of input data (OTUs from Day -1, 0, and 1) and make plots and supplemental table of the top 20 OTUs that were contributing to the 3 models.
```
Rscript code/class._60-40_analysis.R
Rscript code/class_interpretation.R
```
Perform statistical analysis of the 3 classification models.
```
Rscript code/compare_models.R
```
Compare classification model OTUs to the OTUs that varied across sources and/or were altered by clindamycin treatment.
```
Rscript code/venn_diagram_comparisons.R
```

Make figures for the paper.
```
Rscript code/figure_1.R
Rscript code/figure_2.R
Rscript code/figure_3.R
Rscript code/figure_4.R
Rscript code/figure_5.R
Rscript code/figure_6.R
Rscript code/figure_7.R
Rscript code/figure_8.R
Rscript code/figure_S1.R
Rscript code/figure_S2.R
Rscript code/figure_S3.R
Rscript code/figure_S4.R
```
Make excel workbook of supplemental tables.
```
Rscript code/tables_supplemental.R
```

Make supplemental movie for the paper. Note conversion of gif generated from code/pcoa_stats.R requires [FFmpeg](http://ffmpeg.org/). Install FFmpeg using homebrew.
```
brew install ffmpeg
```
Convert gif into .mov file and simultaneously move and rename file into submission directory.
```
ffmpeg -i results/pcoa_over_time.gif -movflags faststart -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" results/pcoa_over_time.mov
cp results/pcoa_over_time.mov submission/movie_S1.mov
```

#### Generate the paper.
```
open submission/manuscript.Rmd and knit to Word or PDF.
```
