AppendList2CSv <- function(l,csvFile) {
  #--------------------------------------------------------------------
  # Write a LIST to CSV file
  #--------------------------------------------------------------------
  out_file <- file(csvFile, open="a")  #creates a file in append mode 
  for (i in seq_along(l)){ 
    # writes the name of the list elements ("A", "B", etc.)
    write.table(names(l)[i],file=out_file,sep=",",dec=".",quote=F,col.names=F,row.names=F) 
    write.table(l[[i]],     file=out_file,sep=",",dec=".",quote=F,col.names=NA,row.names=T) #writes the data.frames 
  } 
  close(out_file)  #close connection to file.csv 
}  
#----------------------------------------------------------------------------------------------------------------------

AppendList2txt <- function(l,csvFile) {
  #--------------------------------------------------------------------
  # Write a LIST to TXT file
  #--------------------------------------------------------------------
  out_file <- file(csvFile, open="a")  #creates a file in append mode 
  for (i in seq_along(l)){ 
    #writes the name of the list elements ("A", "B", etc) 
    write.table(names(l)[i],file=out_file,sep=" ",dec=".",quote=F,col.names=F, row.names=F)
    write.table(l[[i]],     file=out_file,sep=" ",dec=".",quote=F,col.names=NA,row.names=T) #writes the data.frames 
  } 
  close(out_file)  #close connection to file.csv 
}  

DsSplit <- function(ds,trainFrac=3/4,fDet=FALSE,PathDataSet="",iSeed) {
  # ===============================================
  # Dataset spliting in Training and Test (Step 6)
  # ===============================================
  # Inputs
  # - ds = frame dataset object
  # - fDet = flag for detais (TRUE/FALSE)
  # - PathDataSet = pathway for results
  # Output = training and test datasets (to be used for regressions in other functions)
  
  # if datails = TRUE, output files will be created
  my.datf<- ds
  
  # create TRAIN and TEST sets to build a model
  set.seed(iSeed)
  inTrain <- createDataPartition(1:dim(my.datf)[1],p = trainFrac,list = FALSE,groups=2)
  # groups==2 forces to NOT partition
  # based on quantiles of numeric values
  
  my.datf.train<- my.datf[inTrain,]            # TRAIN dataset frame         
  my.datf.test <- my.datf[-inTrain,]           # TEST dataset frame
  
  if (fDet == TRUE) {
    # write the TRAIN and TEST set files
    # the index of each row will in the dataset will not be saved (row.names=F)
    outTrain <- file.path(PathDataSet,paste("ds.Train.split",iSeed,".csv")) # the same folder as the input
    write.csv(my.datf.train,outTrain,row.names=FALSE)
    outTest <- file.path(PathDataSet,paste("ds.Test.split",iSeed,".csv")) # the same folder as the input
    write.csv(my.datf.test,outTest,row.names=FALSE) 
  }
  MyList<- list("train"=my.datf.train, "test"=my.datf.test) 
  return(MyList)  # return a list with training and test datasets
}

Yrandom<- function(dss,trainFrac,best.reg,best.R2.ts,noYrand,ResBestF,rfe_SVM_param_c,rfe_SVM_param_eps){
  #================================================
  # Y-randomization for the best model (Step 12)
  #================================================
  #    - 1 splitting, 1 CV type, best method
  #    - best.R2.ts will be compared with Yrand.R2.ts
  #    - returns ratios DiffsR2/bestR2
  #   (ex param for rbfDDa: negThrStep)
  # --------------------------------------------------
  
  cat("-> Best model Y-Randomization ...\n")
  dss[,1] <- sample(dss[,1]) # randomize Y values for the entire dataset
  
  # splitting dataset in training and test
  #---------------------------------------
  Yrand.R2.ts <- NULL     # all values of R2 for each Y randomization
  for (i in 1:noYrand){
    iSeed=i               
    dsList  <- DsSplit(dss,trainFrac,F,PathDataSet,iSeed) # return a list with 2 datasets = dsList$train, dsList$test
    # get train and test from the resulted list
    ds.train<- dsList$train
    ds.test <- dsList$test
    
    # Run the caret function with the method from the best method
    #    for one training-test split only; no details, we need only R2 values
    if (best.reg=="lm") {
      my.stats.reg  <- LMreg(ds.train,ds.test,"repeatedcv",i,F,ResBestF)$stat.values # run GLM for each CV and regr method
    }
    if (best.reg=="glmStepAIC") {
      my.stats.reg  <- GLMreg(ds.train,ds.test,"repeatedcv",i,F,ResBestF)$stat.values # run GLM for each CV and regr method
    }
    if (best.reg=="pls") {
      my.stats.reg  <- PLSreg(ds.train,ds.test,"repeatedcv",i,F,ResBestF)$stat.values # run SVRM Radial for each CV and regr method
    }
    if (best.reg=="lasso.RMSE") {
      my.stats.reg  <- LASSOreg(ds.train,ds.test,"repeatedcv",i,F,ResBestF)$stat.values # run SVRM Radial for each CV and regr method
    }
    if (best.reg=="svmRadial") {  
      my.stats.reg  <- SVRMreg(ds.train,ds.test,"repeatedcv",i,F,ResBestF,rfe_SVM_param_c)$stat.values # run SVRM Radial for each CV and regr method
    }
    if (best.reg=="nnet") {  
      my.stats.reg  <- NNreg(ds.train,ds.test,"repeatedcv",i,F,ResBestF)$stat.values # run NNet for each CV and regr method
    } 
    if (best.reg=="rf") {  
      my.stats.reg  <- RFreg(ds.train,ds.test,"repeatedcv",i,F,ResBestF)$stat.values # run NNet for each CV and regr method
    } 
    if (best.reg=="svmRFE") {  
      my.stats.reg  <- SVMRFEreg(ds.train,ds.test,"repeatedcv",i,F,ResBestF,rfe_SVM_param_c,rfe_SVM_param_eps)$stat.values # run NNet for each CV and regr method
    } 
    if (best.reg=="glmnet") {  
      my.stats.reg  <- ENETreg(ds.train,ds.test,"repeatedcv",i,F,ResBestF)$stat.values # run NNet for each CV and regr method
    }  
    if (best.reg=="rfRFE") {  
      my.stats.reg  <- RFRFEreg(ds.train,ds.test,"repeatedcv",i,F,ResBestF)$stat.values # run NNet for each CV and regr method
    } 
    
    #     if (best.reg=="rbfDDA") {  
    #       my.stats.reg  <- RBF_DDAreg(ds.train,ds.test,"repeatedcv",negThrStep,i,F,ResBestF)$stat.values # run SVRM Radial for each CV and regr method
    #    }
    
    Yrand.R2.ts <- c(Yrand.R2.ts,my.stats.reg$R2.ts) # adding test R2 value Y randomization
  }
  
  R2diffsPerBestR2 <- NULL
  
  if (is.na(my.stats.reg$R2.ts)) { # check for NA values 
    cat("       --> Y-Randomization error due to NA values!\n")
    write("Y-Randomization error due to NA values!", file=ResBestF,append=T)
  }
  else{
    # get histogram for differences between best R2 and the values for each Y randomization
    R2diffs          <- abs(Yrand.R2.ts - best.R2.ts) # absolute differences between R2 values (best model vs Y randomized results)
    R2diffsPerBestR2 <- abs(R2diffs/best.R2.ts)            # the same difference in percents
    
    pdf(file=paste(ResBestF,".Yrand.Hist.pdf",sep=""))    # save histogram if ratio diffs R2 into PDF for Y random
    Yrand.hist  <- hist(R2diffsPerBestR2)         # draw histogram of the ratio diffs/Best R2 for Y random
    dev.off()
    
    write.table("Y randomization test: ",file=ResBestF,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("=====================", file=ResBestF,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Diffs R2 (Best Model - Y rand):",file=ResBestF,append=T,sep=",",col.names=F,row.names=F,quote=F)
    AppendList2CSv(R2diffs, ResBestF)
    write.table("Summary Difs:",file=ResBestF,append=T,sep=",",col.names=F,row.names=F,quote=F)
    AppendList2CSv(summary(R2diffs), ResBestF)
    write.table("Ratio Diffs R2 / Best R2 (Best Model - Y rand):",file=ResBestF,append=T,sep=",",col.names=F,row.names=F,quote=F)
    AppendList2CSv(R2diffsPerBestR2, ResBestF)
    write.table("Summary Difs %:",file=ResBestF,append=T,sep=",",col.names=F,row.names=F,quote=F)
    AppendList2CSv(summary(R2diffsPerBestR2),ResBestF)
  }
  return(R2diffsPerBestR2) # return the ratio of diffs with the best R2 from the same 
}

findResamps.funct<- function(caret.obj){
  #=============================================================================
  # A function to find the number of re-samples for caret, rfe or sbf objects 
  # from caret package 
  #=============================================================================
  #caret.obj== caret object of class train, rfe or sbf 
  
  in.caret.obj<- caret.obj$control$index
  return(length(in.caret.obj))
}

impute.funct<- function(ds,FUN=mean){
  #=============================================================================
  # A function to impute missing values from columns of matrix or data frame
  # using the mean value as the default 
  #=============================================================================
  #ds== data.frame or matrix to be imputed 
  
  sum.na <- apply(ds,2,sum)
  ind.na <- which(is.na(sum.na)!=FALSE)
  
  ds.imputeV <- apply(as.matrix(ds[,ind.na]),2,function(x)FUN(x,na.rm=T))
  ds.imputeI <- apply(as.matrix(ds[,ind.na]),2,function(x)which(is.na(x)))
  
  if(is.list(ds.imputeI)!=TRUE){ds.imputI<- list(ds.imputeI)}
  
  for(i in 1:length(ds.imputI)){ds[ds.imputI[[i]],ind.na[i]]<- ds.imputeV[i]}
  return(ds)
}
