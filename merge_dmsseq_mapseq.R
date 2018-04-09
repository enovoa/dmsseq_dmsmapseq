########################################################
### MERGE DMS-SEQ & DMS-MaPSEQ
########################################################
### Eva Maria Novoa (e.novoa@garvan.org.au)
### Feb 2018 
### Distributed GNU Affero General Public License v3.0
########################################################


mapseq_masking_based_on_prediction<- function(dat) {
	dat$prediction2<-abs(as.numeric(dat$prediction) - 2)
	dat$merge.rep1 = ifelse(dat$prediction2 > 0, 1 ,dat$react.B1)
	dat$merge.rep2 = ifelse(dat$prediction2 > 0, 1 ,dat$react.B2)
	return(dat)	
}

build_roc_plot_smoothed<-function(dat,smooth_choice,label_smooth) {
	library(pROC)	
	dat$paired_status_ROC_dmsseq<-abs(as.numeric(dat$paired_status) - 2)
	roc1<-roc(dat$paired_status_ROC_dmsseq,dat$missmatch_freq.rep1, smooth=smooth_choice) 
	roc2<- roc(dat$paired_status_ROC_dmsseq, dat$missmatch_freq.rep2, smooth=smooth_choice)
	roc3<-roc(dat$paired_status_ROC_dmsseq,dat$react.B1, smooth=smooth_choice) 
	roc4<- roc(dat$paired_status_ROC_dmsseq, dat$react.B2, smooth=smooth_choice)
	roc5<-roc(dat$paired_status_ROC_dmsseq,dat$merge.rep1, smooth=smooth_choice) 
	roc6<- roc(dat$paired_status_ROC_dmsseq, dat$merge.rep2, smooth=smooth_choice)
	
	pdfName=paste("./results/",paste(label_smooth,".ROC.pdf",sep=""),sep="")
	pdf(file=pdfName, height=7, width=7)
	mycols=c("red1","red2","blue1","blue3","green1","green2")
	plot(roc1, col=mycols[1], main="ROC curve",lwd=3)
	plot(roc2,add=TRUE, col=mycols[2],lwd=3)
	plot(roc3,add=TRUE, col=mycols[3],lwd=3)
	plot(roc4,add=TRUE, col=mycols[4],lwd=3)
	plot(roc5,add=TRUE, col=mycols[5],lwd=3)
	plot(roc6,add=TRUE, col=mycols[6],lwd=3)
	legend('bottomright', cex=0.7,ncol=1,c("MaP.rep1","MaP.rep2","RT_stop.rep1","RT_stop.rep2","merged RTstop+MaP.rep1","merged RTstop+MaP.rep2"),col=mycols,lwd=3, bty = "n")
	dev.off()
}

run_merge_dmsseq_mapseq<- function(input) {
	dat<-read.table(input)
	dat.masked<-mapseq_masking_based_on_prediction(dat)
	print("Writing output file: ./results/mismatches_and_accessibility_merged.txt")
	write.table(dat.masked,file="./results/mismatches_and_accessibility_merged.txt",sep="\t",quote=F)
	print("Writing ROC curve plot: ./results/smoothed.ROC.pdf")
	build_roc_plot_smoothed(dat.masked,TRUE,"smoothed")
}

#TO RUN
#run_merge_dmsseq_mapseq("sample_data/mismatches_and_accessibility_for_merging.txt")
