########################################################
### FILTER MaP
########################################################
### Eva Maria Novoa (e.novoa@garvan.org.au)
### Feb 2018 
### Distributed GNU Affero General Public License v3.0
########################################################

# Read mismatch file
from_file_to_clean_dat<-function(myfile) {

	dat<-read.table(myfile) 
	colnames(dat)<-c("chr","pos","strand","ref_nuc","A","C","G","T","num_miss")
	dat$coverage<-rowSums(dat[,5:8])
	dat$missmatch_freq<-dat$num_miss/dat$coverage
	dat$pos<-dat$pos+1
	row.names(dat)<-paste(dat$chr,dat$pos,sep=".")
	return(dat)
}

# Filtering functions
from_alldat_to_filtered<- function(dat) {
	# 1. Filter based on % missmatch and min_num_missmatch
	dat<-dat[dat$coverage >= 50,] 		 #  Minimum coverage of 50
	dat_missmatched<-dat[dat$missmatch_freq < 0.25,] #  Filter for missmatch_frq
	dat_missmatched<-dat_missmatched[dat_missmatched$num_miss >=2 ,]		 #  Minimum of 2 SNPs
	return (dat_missmatched)
}

reshape_merged_mismatch_datasets<-function(dat1, label1, label2) {
	dat<-dat1[,2:dim(dat1)[2]]
	row.names(dat)<-dat1$Row.names
	colnames(dat)<-paste(rep(colnames(dat1),2),c(rep(label1,11),rep(label2,11)),sep=".")
	print (head(dat))
	return(dat)
}



## To execute:

## 1. Read data
# mismatches_dat.rep1<-from_file_to_clean_dat("mismatches.rep1.txt")  
# mismatches_dat.rep2<-from_file_to_clean_dat("mismatches.rep2.txt")   
## 2. Filter by coverage and mismatch frequency 
# mismatches_dat.rep1.no_snps<-from_alldat_to_filtered(mismatches_dat.rep1)
# mismatches_dat.rep2.no_snps<-from_alldat_to_filtered(mismatches_dat.rep2)
## 3. Filter by replicable mismatches
# mismatches_dat.no_snps.reps<-merge(mismatches_dat.rep1.no_snps,mismatches_dat.rep2.no_snp, by=0, all=TRUE))
# mismatches_dat.no_snps.reps.CLEAN<-reshape_merged_mismatch_datasets(mismatches_dat.no_snps.reps,"rep1","rep2")
## 4. Filter RNASeq positions

