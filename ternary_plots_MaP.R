########################################################
### BUILD TERNARY PLOTS
########################################################
### Eva Maria Novoa (e.novoa@garvan.org.au)
### Feb 2018
### Distributed GNU Affero General Public License v3.0
########################################################

library("ggtern")

myggtern_plot_A<-function(mydf,C,G,T,mytitle) {
	plot<-ggtern(mydf,aes(C,G,T))+ geom_point(aes(color = pair_all.rep2), mydf) + ggtitle(mytitle)
	return(plot)
	#plot+ theme_bw() + theme_nogrid()
}

myggtern_plot_C<-function(mydf,A,G,T,mytitle) {
	plot<-ggtern(mydf,aes(A,G,T))+ geom_point(aes(color = pair_all.rep2), mydf) + ggtitle(mytitle)
	return(plot)
	#plot+ theme_bw() + theme_nogrid()
}

myggtern_plot_for_As_Cs<-function(dat_freq,sampleName){
	dat.A<-dat_freq[dat_freq$ref_nuc.rep1=="A",]
	dat.C<-dat_freq[dat_freq$ref_nuc.rep1=="C",]

	pdf(paste(sampleName,"TERN_PLOT.pdf",sep="."),height=7,width=7)
	plot(myggtern_plot_A(dat.A,C,G,T,paste("A missmatches",sampleName,sep=", ")))
	plot(myggtern_plot_C(dat.C,A,G,T,paste("C missmatches",sampleName,sep=", ")))
	dev.off()
}

run_ternary_plots<-function(dat_freq,dat_freqName) {
	dat_freq$A<-round(rowMeans(dat_freq[,c("A.rep1","A.rep2")]))
	dat_freq$C<-round(rowMeans(dat_freq[,c("C.rep1","C.rep2")]))
	dat_freq$G<-round(rowMeans(dat_freq[,c("G.rep1","G.rep2")]))
	dat_freq$T<-round(rowMeans(dat_freq[,c("T.rep1","T.rep2")]))
	myggtern_plot_for_As_Cs(na.omit(dat_freq),dat_freqName)
}


## To execute:
# run_ternary_plots(input_file,"title")
# Sample input_file: 
