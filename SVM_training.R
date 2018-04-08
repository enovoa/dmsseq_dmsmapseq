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
	#newdat$prediction2<-abs(as.numeric(newdat$prediction) - 2)
	
	print(summary(svm_model))
	print(table(newdat$prediction,newdat$paired_status))

	return(svm_model)
}

test_svm_return_Pred<-function(dat,svm_model,selected_features) {
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
	
	#print(head(newdat))
	#newdat$prediction2<-abs(as.numeric(newdat$prediction) - 2)
	
	print("---> This is the summary of the SVM model")
	print(summary(svm_model))
	print("---> This is the contingency table of the experimental pairing status versus predicted pairing status")
	print(table(newdat$prediction,newdat$paired_status))

	return(newdat)
}

run_SVM_training<-function(dat_input, output_file) {
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
    dat.withPred<-test_svm_return_Pred(dat.TRAIN.As, svm_model.TRAIN.As,c("missmatch_freq.rep1","TGratio.rep1","missmatch_freq.rep2","TGratio.rep2","react.B1","react.B2"))

  	# 6. Subset data and write table
 	write.table(dat.withPred, file="results/mismatches_with_SVM_predictions.txt", quote=F)

}

## SCRIPT:
# run_SVM_training("sample_data/mismatches_for_SVM_training.txt")


