########################################################
### FILTER MaP
########################################################
### Eva Maria Novoa (e.novoa@garvan.org.au)
### Feb 2018 
### Distributed GNU Affero General Public License v3.0
########################################################

## FUNCTIONS:
from_file_to_clean_dat<-function(myfile) {

	dat<-read.table(myfile) 
	colnames(dat)<-c("chr","pos","strand","ref_nuc","A","C","G","T","num_miss")
	dat$coverage<-rowSums(dat[,5:8])
	dat$missmatch_freq<-dat$num_miss/dat$coverage
	dat$pos<-dat$pos+1
	row.names(dat)<-paste(dat$chr,dat$pos,sep=".")
	return(dat)
}

from_alldat_to_filtered<- function(dat) {
	dat<-dat[dat$coverage >= 50,] 		 #  Minimum coverage of 50
	dat_missmatched<-dat[dat$missmatch_freq < 0.25,] #  Filter for missmatch_frq
	dat_missmatched<-dat_missmatched[dat_missmatched$num_miss >=2 ,]		 #  Minimum of 2 SNPs
	return (dat_missmatched)
}

reshape_merged_mismatch_datasets<-function(dat1,dat2, label1, label2) {
	dat<-dat1[,2:dim(dat1)[2]]
	row.names(dat)<-dat1$Row.names
	colnames(dat)<-paste(rep(colnames(dat2),2),c(rep(label1,11),rep(label2,11)),sep=".")
	print (head(dat))
	return(dat)
}



## SCRIPT:


## 1. Read data
mismatches_dmsseq.rep1<-from_file_to_clean_dat("sample_data/mismatches_dmsseq.rep1.txt")
mismatches_dmsseq.rep2<-from_file_to_clean_dat("sample_data/mismatches_dmsseq.rep2.txt")
mismatches_rnaseq.rep1<-from_file_to_clean_dat("sample_data/mismatches_rnaseq.rep1.txt")
mismatches_rnaseq.rep2<-from_file_to_clean_dat("sample_data/mismatches_rnaseq.rep2.txt")

## 2. Filter by coverage and mismatch frequency 
mismatches_dmsseq.rep1.no_snps<-from_alldat_to_filtered(mismatches_dmsseq.rep1)
mismatches_dmsseq.rep2.no_snps<-from_alldat_to_filtered(mismatches_dmsseq.rep2)
mismatches_rnaseq.rep1.no_snps<-from_alldat_to_filtered(mismatches_rnaseq.rep1)
mismatches_rnaseq.rep2.no_snps<-from_alldat_to_filtered(mismatches_rnaseq.rep2)

## 3. Filter by replicable mismatches
mismatches_dmsseq.no_snps.reps<-merge(mismatches_dmsseq.rep1.no_snps,mismatches_dmsseq.rep2.no_snps, by=0, all=FALSE)
mismatches_rnaseq.no_snps.reps<-merge(mismatches_rnaseq.rep1.no_snps,mismatches_rnaseq.rep2.no_snps, by=0, all=FALSE)

## 4. Reshaping, ordering
mismatches_dmsseq.no_snps.reps.CLEAN<-reshape_merged_mismatch_datasets(mismatches_dmsseq.no_snps.reps,mismatches_dmsseq.rep1.no_snps,"rep1","rep2")
mismatches_rnaseq.no_snps.reps.CLEAN<-reshape_merged_mismatch_datasets(mismatches_rnaseq.no_snps.reps,mismatches_rnaseq.rep1.no_snps,"rep1","rep2")

mismatches_dmsseq.no_snps.reps.CLEAN<-mismatches_dmsseq.no_snps.reps.CLEAN[order(mismatches_dmsseq.no_snps.reps.CLEAN$pos.rep1),]
mismatches_rnaseq.no_snps.reps.CLEAN<-mismatches_rnaseq.no_snps.reps.CLEAN[order(mismatches_rnaseq.no_snps.reps.CLEAN$pos.rep1),]

## 5. Filter RNASeq positions
mismatches_dmsseq.rnafilter<-mismatches_dmsseq.no_snps.reps.CLEAN[ ! row.names(mismatches_dmsseq.no_snps.reps.CLEAN) %in% row.names(mismatches_rnaseq.no_snps.reps.CLEAN), ]

## 6. Write table
write.table(mismatches_dmsseq.rnafilter, file="results/mismatches_dmsseq.filtered.txt", quote=F)


