########################################################
### BUILD TERNARY PLOTS
########################################################
### Eva Maria Novoa (e.novoa@garvan.org.au)
### Feb 2018
### Distributed GNU Affero General Public License v3.0
########################################################



myggtern_plot_A<-function(mydf,C,G,T,mytitle) {
    library("ggtern")
	plot<-ggtern(mydf,aes(C,G,T))+ geom_point(aes(color = log(mean_coverage)), mydf) + ggtitle(mytitle)
	return(plot)
}

myggtern_plot_C<-function(mydf,A,G,T,mytitle) {
    library("ggtern")
	plot<-ggtern(mydf,aes(A,G,T))+ geom_point(aes(color = log(mean_coverage)), mydf) + ggtitle(mytitle)
	return(plot)
}

myggtern_plot_for_As_Cs<-function(dat_freq,sampleName){
    library("ggtern")
	dat.A<-dat_freq[dat_freq$ref_nuc.rep1=="A",]
	dat.C<-dat_freq[dat_freq$ref_nuc.rep1=="C",]

	pdf(paste("./results",paste(sampleName,"TERN_PLOT.pdf",sep="."),sep="/"),height=7,width=7)
	plot(myggtern_plot_A(dat.A,C,G,T,paste("A mismatches",sampleName,sep=", ")))
	plot(myggtern_plot_C(dat.C,A,G,T,paste("C mismatches",sampleName,sep=", ")))
	dev.off()
}

run_plot_4_nucleotides<-function(dat_freq,dat_freqName) {
    library("ggtern")
	dat_freq$A<-round(rowMeans(dat_freq[,c("A.rep1","A.rep2")]))
	dat_freq$C<-round(rowMeans(dat_freq[,c("C.rep1","C.rep2")]))
	dat_freq$G<-round(rowMeans(dat_freq[,c("G.rep1","G.rep2")]))
	dat_freq$T<-round(rowMeans(dat_freq[,c("T.rep1","T.rep2")]))
	myggtern_plot_for_As_Cs(na.omit(dat_freq),dat_freqName)
}

run_ternary_plot<-function(dat_input) {
    dat<-read.table(dat_input)
    run_plot_4_nucleotides(dat,"ternary_plot")
}

# SCRIPT:
# run_ternary_plot("sample_data/mismatches_for_ternary_plot.txt")
