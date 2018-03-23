# dmsseq_dmsmapseq
Analysis of DMSSeq and DMSMaPSeq datasets

Here we provide several scripts that can facilitate the filtering, processing and merging of DMS-Seq and DMS-MaPSeq datasets. 

A small dataset to demo the scripts is included (folder: sample_data)
Expected output of the scripts is also included (folder: results)

# 1. System requirements and installation
These scripts are written in R, and have been tested under R-3.3.2. 
R versions should be at least XXX due to specific R package dependency.

Packages that need to be installed:
- ggtern (version: 2.2.1 or greater)
- caTools (version 1.17.1 or greater)
- e1071 (version 1.6-8)

To install the R packages, start R, and prompt:
install.packages("ggtern")
install.packages("caTools")
install.packages("e1071")

No prior installation of additional software is required.

# 2. Instructions for use

* To run filtering of MaP data: 
> source("filter_MaP.R")
> run_MaP_filtering("sample_data/mismatches_dmsseq.rep1.txt","sample_data/mismatches_dmsseq.rep2.txt","sample_data/mismatches_rnaseq.rep1.txt","sample_data/mismatches_rnaseq.rep2.txt","results/mismatches_dmsseq.filtered.txt")

* To build ternary plots: 
> source("ternary_plots_MaP.R")
> run_ternary_plot("sample_data/mismatches_for_ternary_plot.txt")

* To build SVM and get predictions:
> source("SVM_training.R")
> run_SVM_training("sample_data/mismatches_for_SVM_training.txt")

* To merge DMS-Seq (RT-drop-off) and MaPSeq (mutational profiling) data:
> source("merge_dmsseq_mapseq.R")
> run_merge_dmsseq_mapseq("sample_data/mismatches_and_accessibility_for_merging.txt")



