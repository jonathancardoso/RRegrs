RFreg <- function(my.datf.train,my.datf.test,sCV,iSplit=1,fDet=F,outFile="") {
  #======================================
  # Basic RandomForest
  #======================================
  
  net.c = my.datf.train[,1]   # make available the names of variables from training dataset
  RegrMethod <- "rf" # type of regression
  
  # Define the CV conditions
  ctrl<- trainControl(method=sCV, number=10,repeats=10,#number=10,repeats=10,
                      summaryFunction=defaultSummary)
  
  tuneParam = data.frame(.mtry=c(ncol(my.datf.train)/3,ncol(my.datf.train)/2,ncol(my.datf.train)))
  # Train the model using only training set
  set.seed(iSplit)
  rf.fit<- train(net.c~.,data=my.datf.train,
                 method='rf', trControl=ctrl,
                 metric='RMSE',ntree=1500,tuneGrid =tuneParam)
  
  #------------------------------
  # Training RESULTS
  #------------------------------
  RMSE.tr  <- rf.fit$results[rownames(rf.fit$bestTune),2]
  R2.tr    <- rf.fit$results[rownames(rf.fit$bestTune),3]
  if (sCV == "repeatedcv"){ # if 10-fold CV
    RMSEsd.tr <- rf.fit$results[rownames(rf.fit$bestTune),4]
    R2sd.tr   <- rf.fit$results[rownames(rf.fit$bestTune),5]
  }
  if (sCV == "LOOCV"){ # if LOOCV
    RMSEsd.tr <- 0 # formulas will be added later! 
    R2sd.tr   <- 0 # formulas will be added later!
  }
  
  #------------------------------------------------
  # RMSE & R^2, for train/test respectively
  #------------------------------------------------
  rf.train.res <- getTrainPerf(rf.fit)
  rf.test.res  <- postResample(predict(rf.fit,my.datf.test),my.datf.test[,1])
  
  #------------------------------------------------
  # Adj R2, Pearson correlation
  #------------------------------------------------
  pred.tr     <- predict(rf.fit,my.datf.train) # predicted Y for training
  pred.ts     <- predict(rf.fit,my.datf.test)  # predicted Y for test
  noFeats.fit <- length(predictors(rf.fit))    # no. of features from the fitted model
  Feats.fit   <- paste(predictors(rf.fit),collapse="+") # string with the features included in the fitted model
  
  ds.full     <- rbind(my.datf.train,my.datf.test)
  pred.both   <- predict(rf.fit,ds.full)       # predicted Y
  mae.tr      <- mae(my.datf.train[,1],pred.tr)
  mae.ts      <- mae(my.datf.test[,1],pred.ts)
  adjR2.tr    <- r2.adj.funct(my.datf.train[,1],pred.tr,noFeats.fit)
  adjR2.ts    <- r2.adj.funct(my.datf.test[,1],pred.ts,noFeats.fit)
  corP.ts     <- cor(my.datf.test[,1],pred.ts)
  
  adjR2.both  <- r2.adj.funct(ds.full[,1],pred.both,noFeats.fit)
  RMSE.both   <- rmse.funct(ds.full[,1],pred.both)
  r2.both     <- r2.funct(ds.full[,1],pred.both)
  
  # Generate the output list with statistics for each cross-validation type
  # ------------------------------------------------------------------------------
  my.stats <- list("RegrMeth"     = RegrMethod,
                   "Split No"     = as.numeric(iSplit),     # from function param
                   "CVtype"       = sCV,                    # from function param
                   "NoModelFeats" = as.numeric(noFeats.fit),
                   "ModelFeats"   = Feats.fit,
                   "adjR2.tr"  = as.numeric(adjR2.tr),
                   "RMSE.tr"   = as.numeric(RMSE.tr),
                   "MAE.tr"    = as.numeric(mae.tr),
                   "R2.tr"     = as.numeric(R2.tr),
                   "RMSEsd.tr" = as.numeric(RMSEsd.tr),
                   "R2sd.tr"   = as.numeric(R2sd.tr),
                   "adjR2.ts"= as.numeric(adjR2.ts),
                   "RMSE.ts" = as.numeric((rf.test.res["RMSE"])),
                   "MAE.ts"    = as.numeric(mae.ts),
                   "R2.ts"   = as.numeric((rf.test.res["Rsquared"])),
                   "corP.ts" = as.numeric(corP.ts),
                   "adjR2.both" = as.numeric(adjR2.both),
                   "RMSE.both"  = as.numeric(RMSE.both),
                   "R2.both"    = as.numeric(r2.both))
  
  #---------------------------------------------------------------------
  # Write to file DETAILS for GLM for each cross-validation method
  #---------------------------------------------------------------------
  if (fDet==T) {   # if flag for details if T, print details about any resut
    write("RRegr package | eNanoMapper", file=outFile, append=T)
    write.table(paste("Regression method: ", RegrMethod), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("Split no.: ", iSplit), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("CV type: ", sCV),      file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Training Set Summary: ",     file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(summary(my.datf.train),       file=outFile,append=T,sep=",",col.names=T,quote=F)
    write.table("Test Set Summary: ",         file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(summary(my.datf.test),        file=outFile,append=T,sep=",",col.names=T,quote=F)   
    
    
    write.table("Predictors: ",     file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(rf.fit), file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Trainig Results: ",     file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(rf.train.res),file=outFile,append=T,sep=",",col.names=T,quote=F)
    write.table("Test Results: ",        file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(rf.test.res), file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Full Statistics: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(my.stats,            file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    # Variable Importance (max top 20)
    FeatImp <- importance(rf.fit$finalModel, scale = T)
    FeatImp = FeatImp[order(FeatImp,decreasing=T),]
    components = length(FeatImp)  # default plot all feature importance
    if (length(FeatImp)>20){     # if the number of features is greater than 20, use only 20
      components = 20
    }
    # Append feature importance to output details
    AppendList2CSv(FeatImp,outFile)
    
    fitModel <- rf.fit$finalModel
    
    # =============================================================================
    # Assessment of Applicability Domain (plot leverage)
    # =============================================================================
    
    # Residuals
    resids <- fitModel$predicted-fitModel$y # residuals
    write.table("Residuals of the fitted model: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(data.frame(resids), file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write residuals
    
    predVals.pls.ad <- pred.ts
    Traind.pls= as.matrix(my.datf.train)
    Testd.pls = as.matrix(my.datf.test)
    mat.Traind.pls<- t(Traind.pls) %*%(Traind.pls) 
    det.Traind.pls<- det(mat.Traind.pls)
    
    if(det.Traind.pls!=0){
      Hat.train = diag(Traind.pls %*% solve(t(Traind.pls) %*%(Traind.pls), tol=1e-40)  %*% t(Traind.pls))
      Hat.test  = diag(Testd.pls  %*% solve(t(Traind.pls) %*%(Traind.pls), tol=1e-40)  %*% t(Testd.pls))  
      
      # Leverage / Hat values
      hat.fit <- Hat.test          # hat values
      hat.fit.df <- as.data.frame(hat.fit)    # hat data frame
      hat.mean <- mean(hat.fit)               # mean hat values
      hat.fit.df$warn <- ifelse(hat.fit.df[, 'hat.fit']>3*hat.mean, 'x3',ifelse(hat.fit.df[, 'hat.fit']>2*hat.mean, 'x2', '-' ))
      
      write.table("Leverage output: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
      write.table(paste("Mean of hat values: ", hat.mean), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
      write.table("Leverage / Hat values with warnings (X3 & X2 = values 3 & 2 times than hat mean): ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
      write.table(hat.fit.df, file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write hat values and the levels X3, X2 (of hat mean)
      
      #THRESHOLD values: 3m/n, where m is the number of parameters, and n number of observations
      thresh.lever<- (3*(dim(my.datf.train)[2]-1))/dim(my.datf.train)[1] # leverage thresh
      hat.problems<- data.frame(hat.fit[hat.fit>thresh.lever]) # points with high leverage
      
      write.table(paste("Leverage Threshold: ", thresh.lever), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
      write.table("Points with leverage > threshold: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
      write.table(hat.problems, file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F)
    } 
    # PDF with 12 plots
    # --------------------------------------------------------------
    pdf(file=paste(outFile,".",sCV,".","split",iSplit,".pdf",sep=""))
    # par(mfrow = c(3, 4)) # all plots into one page!
    
    plot(my.datf.train[,1],pred.tr,xlab="Yobs", ylab="Ypred", type="b", main="Train Yobs-Ypred") # plot 1
    plot(my.datf.test[,1], pred.ts,xlab="Yobs", ylab="Ypred", type="b", main="Test Yobs-Ypred")  # plot 2
    dotchart(as.matrix(FeatImp),main="Feature Importance")                            # plot 3
    
    # Fitted vs Residuals - plot 4
    plot(fitModel$predicted,resids,
         main="Fitted vs. Residuals for Fitted Model",
         xlab="Fitted", ylab="Residuals")
    abline(h = 0, lty = 2)
    
    # Leverage plots - plot 5
    if(det.Traind.pls!=0){
      plot(hat.fit, type = "h",
           main="Leverage for Fitted Model",
           xlab="Index", ylab="Hat")
      abline(h = thresh.lever, lty = 2, col="red") # leverage thresh
    }
    # plot(FeatImp, top = components,main="Feature Importance") # ERROR !
    dev.off()
    # --------------------------------------------------------------
  }
  
  return(list(stat.values=my.stats, model=rf.fit))  # return a list with statistics and the full model
}
#----------------------------------------------------------------------------------------------------------------------

RFRFEreg <- function(my.datf.train,my.datf.test,sCV,iSplit=1,fDet=F,outFile="") {
  #=============================
  # Random Forest-RFE
  #=============================
  
  net.c = my.datf.train[,1]   # make available the names of variables from training dataset
  RegrMethod <- "rfRFE" # type of regression
  
  # Define the CV conditions
  ctrl<- trainControl(method=sCV, number=5,repeats=1,#number=10,repeats=10,
                      summaryFunction=defaultSummary,verboseIter = F)
  
  rfeCtr = rfeControl(functions = rfFuncs,method="cv",number=10,repeats=10, saveDetails = T, verbose=T,rerank = F,allowParallel=T)#number=10,repeats=10,
  
  sizes = 2^(1:sqrt(ncol(my.datf.train)-1))
  
  # Train the model using only training set
  set.seed(iSplit)
  
  input = as.matrix(my.datf.train[,2:ncol(my.datf.train)])
  rferf.fit = rfe(input,net.c,sizes = sizes,rfeControl=rfeCtr,prob.model =F,trControl=ctrl ,allowParallel=T, tuneGrid=expand.grid(.mtry=c(floor(sqrt(ncol(input))),ncol(input))), metric='RMSE')
  
  #------------------------------
  # Training RESULTS
  #------------------------------
  RMSE.tr  <- rferf.fit$results[rferf.fit$results$Variables ==rferf.fit$bestSubset,2]
  R2.tr    <- rferf.fit$results[rferf.fit$results$Variables ==rferf.fit$bestSubset,3]
  if (sCV == "repeatedcv"){ # if 10-fold CV
    RMSEsd.tr <- rferf.fit$results[rferf.fit$results$Variables ==rferf.fit$bestSubset,4]
    R2sd.tr   <- rferf.fit$results[rferf.fit$results$Variables ==rferf.fit$bestSubset,5]
  }
  if (sCV == "LOOCV"){ # if LOOCV
    RMSEsd.tr <- 0 # formulas will be added later!  TODOOOOOOOOOOOOOO
    R2sd.tr   <- 0 # formulas will be added later!
  }
  
  #------------------------------------------------
  # RMSE & R^2, for train/test respectively
  #------------------------------------------------
  rfesvm.train.res <- rferf.fit$results[ rferf.fit$results$Variables== rferf.fit$bestSubset, c(2,3)]
  rfesvm.test.res  <- postResample(predict(rferf.fit,my.datf.test),my.datf.test[,1])
  
  #------------------------------------------------
  # Adj R2, Pearson correlation
  #------------------------------------------------
  pred.tr     <- predict(rferf.fit,my.datf.train) # predicted Y for training
  pred.ts     <- predict(rferf.fit,my.datf.test)  # predicted Y for test
  noFeats.fit <- length(predictors(rferf.fit))    # no. of features from the fitted model
  Feats.fit   <- paste(predictors(rferf.fit),collapse="+") # string with the features included in the fitted model
  
  ds.full     <- rbind(my.datf.train,my.datf.test)
  pred.both   <- predict(rferf.fit,ds.full)       # predicted Y
  mae.tr      <- mae(my.datf.train[,1],pred.tr)
  mae.ts      <- mae(my.datf.test[,1],pred.ts)
  adjR2.tr    <- r2.adj.funct(my.datf.train[,1],pred.tr,noFeats.fit)
  adjR2.ts    <- r2.adj.funct(my.datf.test[,1],pred.ts,noFeats.fit)
  corP.ts     <- cor(my.datf.test[,1],pred.ts)
  
  adjR2.both  <- r2.adj.funct(ds.full[,1],pred.both,noFeats.fit)
  RMSE.both   <- rmse.funct(ds.full[,1],pred.both)
  r2.both     <- r2.funct(ds.full[,1],pred.both)
  
  # Generate the output list with statistics for each cross-validation type
  # ----------------------------------------------------------------------------
  my.stats <- list("RegrMeth"     = RegrMethod,
                   "Split No"     = as.numeric(iSplit),     # from function param
                   "CVtype"       = sCV,                    # from function param
                   "NoModelFeats" = as.numeric(noFeats.fit),
                   "ModelFeats"   = Feats.fit,
                   "adjR2.tr"  = as.numeric(adjR2.tr),
                   "RMSE.tr"   = as.numeric(RMSE.tr),
                   "MAE.tr"    = as.numeric(mae.tr),
                   "R2.tr"     = as.numeric(R2.tr),
                   "RMSEsd.tr" = as.numeric(RMSEsd.tr),
                   "R2sd.tr"   = as.numeric(R2sd.tr),
                   "adjR2.ts"= as.numeric(adjR2.ts),
                   "RMSE.ts" = as.numeric((rfesvm.test.res["RMSE"])),
                   "MAE.ts"    = as.numeric(mae.ts),
                   "R2.ts"   = as.numeric((rfesvm.test.res["Rsquared"])),
                   "corP.ts" = as.numeric(corP.ts),
                   "adjR2.both" = as.numeric(adjR2.both),
                   "RMSE.both"  = as.numeric(RMSE.both),
                   "R2.both"    = as.numeric(r2.both))
  
  #---------------------------------------------------------------------
  # Write to file DETAILS for GLM for each cross-validation method
  #---------------------------------------------------------------------
  if (fDet==T) {   # if flag for details if T, print details about any resut
    write("RRegr package | eNanoMapper", file=outFile, append=T)
    write.table(paste("Regression method: ", RegrMethod), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("Split no.: ", iSplit), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("CV type: ", sCV),      file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Training Set Summary: ",     file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(summary(my.datf.train),       file=outFile,append=T,sep=",",col.names=T,quote=F)
    write.table("Test Set Summary: ",         file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(summary(my.datf.test),        file=outFile,append=T,sep=",",col.names=T,quote=F)   
    
    
    write.table("Predictors: ",     file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(rferf.fit), file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Trainig Results: ",     file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(rfesvm.train.res,file=outFile,append=T,sep=",",col.names=T,quote=F)
    write.table("Test Results: ",        file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(rfesvm.test.res, file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Full Statistics: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(my.stats,            file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    # Variable Importance (max top 20)
    FeatImp <- importance(rferf.fit$fit, scale = T)
    FeatImp = FeatImp[order(FeatImp,decreasing=T),]
    
    components = length(FeatImp)  # default plot all feature importance
    if (length(FeatImp)>20){     # if the number of features is greater than 20, use only 20
      components = 20
    }
    # Append feature importance to output details
    AppendList2CSv(FeatImp,outFile)
    
    fitModel <- rferf.fit$fit
    
    # =============================================================================
    # Assessment of Applicability Domain (plot leverage)
    # =============================================================================
    
    # Residuals
    resids <- pred.both-ds.full[,1]    # residuals
    write.table("Residuals of the fitted model: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(data.frame(resids), file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write residuals
    
    predVals.pls.ad <- pred.ts
    Traind.pls= as.matrix(my.datf.train)
    Testd.pls = as.matrix(my.datf.test)
    mat.Traind.pls<- t(Traind.pls) %*%(Traind.pls) 
    det.Traind.pls<- det(mat.Traind.pls)
    
    if(det.Traind.pls!=0){
      
      Hat.train = diag(Traind.pls %*% solve(t(Traind.pls) %*%(Traind.pls), tol=1e-40)  %*% t(Traind.pls))
      Hat.test  = diag(Testd.pls  %*% solve(t(Traind.pls) %*%(Traind.pls), tol=1e-40)  %*% t(Testd.pls))  
      
      # Leverage / Hat values
      hat.fit <- Hat.test          # hat values
      hat.fit.df <- as.data.frame(hat.fit)    # hat data frame
      hat.mean <- mean(hat.fit)               # mean hat values
      hat.fit.df$warn <- ifelse(hat.fit.df[, 'hat.fit']>3*hat.mean, 'x3',ifelse(hat.fit.df[, 'hat.fit']>2*hat.mean, 'x2', '-' ))
      
      write.table("Leverage output: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
      write.table(paste("Mean of hat values: ", hat.mean), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
      write.table("Leverage / Hat values with warnings (X3 & X2 = values 3 & 2 times than hat mean): ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
      write.table(hat.fit.df, file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write hat values and the levels X3, X2 (of hat mean)
      
      #THRESHOLD values: 3m/n, where m is the number of parameters, and n number of observations
      thresh.lever<- (3*(dim(my.datf.train)[2]-1))/dim(my.datf.train)[1] # leverage thresh
      hat.problems<- data.frame(hat.fit[hat.fit>thresh.lever]) # points with high leverage
      
      write.table(paste("Leverage Threshold: ", thresh.lever), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
      write.table("Points with leverage > threshold: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
      write.table(hat.problems, file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F)
    }
    
    # PDF with 12 plots
    # --------------------------------------------------------------
    pdf(file=paste(outFile,".",sCV,".","split",iSplit,".pdf",sep=""))
    # par(mfrow = c(3, 4)) # all plots into one page!
    
    plot(my.datf.train[,1],pred.tr,xlab="Yobs", ylab="Ypred", type="b", main="Train Yobs-Ypred") # plot 1
    plot(my.datf.test[,1], pred.ts,xlab="Yobs", ylab="Ypred", type="b", main="Test Yobs-Ypred")  # plot 2
    
    dotchart(FeatImp,main="Feature Importance") # plot 3
    
    # Fitted vs Residuals - plot 4
    plot(pred.both,resids,
         main="Fitted vs. Residuals for Fitted Model",
         xlab="Fitted", ylab="Residuals")
    abline(h = 0, lty = 2)
    
    # Leverage plots - plot 5
    if(det.Traind.pls!=0){
      plot(hat.fit, type = "h",
           main="Leverage for Fitted Model",
           xlab="Index", ylab="Hat")
      abline(h = thresh.lever, lty = 2, col="red") # leverage thresh
    }
    # plot(FeatImp, top = components,main="Feature Importance") # ERROR !
    dev.off()
    # --------------------------------------------------------------
  }
  
  return(list(stat.values=my.stats, model=rferf.fit))  # return a list with statistics and the full model
}
