########################################################
### SVM training (RNA structure pairing status)
########################################################
### Eva Maria Novoa (e.novoa@garvan.org.au)
### Feb 2018 
### Distributed GNU Affero General Public License v3.0
########################################################


subdivide_training_testing<-function(data, type) {
	# Get reproducible sampling	
	require(caTools)
	set.seed(101) 
	sample = sample.split(data$chr, SplitRatio = 0.25)	
	
	# Test or train
	if (type=="train") {
		train = subset(data, sample == TRUE)
		print (dim(train))
		return(train)
		
	}
	if (type=="test") {
		test  = subset(data, sample == FALSE)		
		print(dim(test))
		return(test)
	}
	
}

train_svm_return_model<-function(dat,selected_features) {
    set.seed(101)
	library("e1071")

	x<-dat[,selected_features]
	y<-as.factor(dat$paired_status)
	svm_model <- svm(x,y, cross=5)
	pred <- predict(svm_model,x)
	newdat<-merge(dat,pred, by=0, all=TRUE)
	row.names(newdat)<-newdat[,1]
	newdat<-newdat[,2:dim(newdat)[2]]
	colnames(newdat)<-c(colnames(dat),"prediction")
	print(head(newdat))
	newdat$prediction2<-abs(as.numeric(newdat$prediction) - 2)
	
	print(summary(svm_model))
	print(table(newdat$prediction,newdat$paired_status))

	return(svm_model)
}

test_svm_return_Pred<-function(dat,selected_features) {
	set.seed(101)
    library("e1071")

	x<-dat[,selected_features]
	y<-as.factor(dat$paired_status)
	svm_model <- svm(x,y, cross=5)
	pred <- predict(svm_model,x)
	newdat<-merge(dat,pred, by=0, all=TRUE)
	row.names(newdat)<-newdat[,1]
	newdat<-newdat[,2:dim(newdat)[2]]
	colnames(newdat)<-c(colnames(dat),"prediction")
	print(head(newdat))
	newdat$prediction2<-abs(as.numeric(newdat$prediction) - 2)
	
	print(summary(svm_model))
	print(table(newdat$prediction,newdat$paired_status))

	return(svm_model)
}

run_SVM_training<-function(dat_input) {
    dat<-read.table(dat_input)

    # 1. Subdivide into training and testing
    dat.TRAIN<-subdivide_training_testing(dat,"train")
    dat.TEST<-subdivide_training_testing(dat,"test")

    # 2. Subdivide into As and Cs
    dat.TRAIN.As<-dat.TRAIN[dat.TRAIN$ref_nuc=="A",]
    dat.TRAIN.Cs<-dat.TRAIN[dat.TRAIN$ref_nuc=="C",]

    # 3. Train SVM
    svm_model.TRAIN.As<-train_svm_return_model(dat.TRAIN.As, c("missmatch_freq.rep1","TGratio.rep1","missmatch_freq.rep2","TGratio.rep2","react.B1","react.B2"))

    # 4. Test SVM
    dat.withPred<-run_svm_return_dataset_withPred(dat.TRAIN.As, svm_model.TRAIN.As,c("missmatch_freq.rep1","TGratio.rep1","missmatch_freq.rep2","TGratio.rep2","react.B1","react.B2"))

    # 5. Rearrange data
    selected<-c("pos.rep1","chr.rep1","strand.rep1","ref_nuc.rep1","A.rep1","C.rep1","G.rep1","T.rep1","num_miss.rep1","coverage.rep1"  ,"missmatch_freq.rep1","A.rep2","C.rep2","G.rep2","T.rep2","num_miss.rep2","coverage.rep2","missmatch_freq.rep2" ,"react.rep1","react.rep2","pair_AC.rep1","TAratio.rep1","TAratio.rep2","TGratio.rep1" ,"TGratio.rep2","react.B1","react.B2")
    newcolnames<-c("pos","chr","strand","ref_nuc","A.rep1","C.rep1","G.rep1","T.rep1","num_miss.rep1","coverage.rep1"  ,"missmatch_freq.rep1","A.rep2","C.rep2","G.rep2","T.rep2","num_miss.rep2","coverage.rep2","missmatch_freq.rep2" ,"react.rep1","react.rep2","paired_status","TAratio.rep1","TAratio.rep2","TGratio.rep1" ,"TGratio.rep2","react.B1","react.B2")
 
    ## FINISH !!! ##
}

## SCRIPT:
# run_SVM_training("sample_data/mismatches_for_SVM_training.txt")


