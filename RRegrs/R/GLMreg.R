GLMreg <- function(my.datf.train,my.datf.test,sCV,iSplit=1,fDet=F,outFile="") {
  #======================================================
  # 8.2- GLM stepwise regression - based on AIC (caret)
  #======================================================
  # Inputs:
  # - my.datf.train,my.datf.test = training and test dataset frames
  # - sCV = type of cross-validation such as repeatedcv, LOOCV, etc.
  # - iSplit = index of splitalse
  # - fDet = flag for detais (True/F)
  # - outFile = output file for GLM details
  # Output:
  # - list of statistics equal with the header introduced in the main script and the full model
  #   (tr = train, ts = test, both = tr+ts = full dataset)
  # -----------------------------------------------------------------------------------------------
  
  #library(caret)
  #attach(my.datf.train)    # make available the names of variables from training dataset
  net.c = my.datf.train[,1] # dependent variable is the first column in Training set
  RegrMethod <- "glmStepAIC" # type of regression
  
  # Define the CV conditions
  ctrl<- trainControl(method=sCV, number=10,repeats=10,
                      summaryFunction=defaultSummary)
  
  # Train the model using only training set
  set.seed(iSplit)
  glm.fit<- train(net.c~.,data=my.datf.train,
                  method='glmStepAIC', tuneLength=10, trControl=ctrl,
                  metric='RMSE')
  
  #------------------------------
  # Training RESULTS
  #------------------------------
  RMSE.tr  <- glm.fit$results[,2]
  R2.tr    <- glm.fit$results[,3]
  if (sCV == "repeatedcv"){ # if 10-fold CV
    RMSEsd.tr <- glm.fit$results[,4]
    R2sd.tr   <- glm.fit$results[,5]
  }
  if (sCV == "LOOCV"){ # if LOOCV
    RMSEsd.tr <- 0 # formulas will be added later!
    R2sd.tr   <- 0 # formulas will be added later!
  }
  
  #------------------------------------------------
  # RMSE & R^2, for train/test respectively
  #------------------------------------------------
  lm.train.res <- if(sCV != "none") { getTrainPerf(glm.fit)}
  lm.test.res  <- postResample(predict(glm.fit,my.datf.test),my.datf.test[,1])
  
  #------------------------------------------------
  # Adj R2, Pearson correlation
  #------------------------------------------------
  pred.tr     <- predict(glm.fit,my.datf.train) # predicted Y
  pred.ts     <- predict(glm.fit,my.datf.test)  # predicted Y
  noFeats.fit <- length(predictors(glm.fit))    # no. of features from the fitted model
  Feats.fit   <- paste(predictors(glm.fit),collapse="+") # string with the features included in the fitted model
  
  ds.full     <- rbind(my.datf.train,my.datf.test)
  pred.both   <- predict(glm.fit,ds.full)       # predicted Y
  mae.tr      <- mae(my.datf.train[,1],pred.tr)
  mae.ts      <- mae(my.datf.test[,1],pred.ts)
  adjR2.tr    <- r2.adj.funct(my.datf.train[,1],pred.tr,noFeats.fit)
  adjR2.ts    <- r2.adj.funct(my.datf.test[,1],pred.ts,noFeats.fit)
  corP.ts     <- cor(my.datf.test[,1],pred.ts)
  
  adjR2.both  <- r2.adj.funct(ds.full[,1],pred.both,noFeats.fit)
  RMSE.both   <- rmse.funct(ds.full[,1],pred.both)
  r2.both     <- r2.funct(ds.full[,1],pred.both)
  
  # Generate the output list with statistics for each cross-validation type
  # -------------------------------------------------------------------------
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
                   "RMSE.ts" = as.numeric((lm.test.res["RMSE"][[1]])),
                   "MAE.ts"    = as.numeric(mae.ts),
                   "R2.ts"   = as.numeric((lm.test.res["Rsquared"][[1]])),
                   "corP.ts" = as.numeric(corP.ts),
                   "adjR2.both" = as.numeric(adjR2.both),
                   "RMSE.both"  = as.numeric(RMSE.both),
                   "R2.both"    = as.numeric(r2.both))
  
  #---------------------------------------------------------------------
  # Write to file DETAILS for GLM for each cross-validation method
  #---------------------------------------------------------------------
  if (fDet==T) {   # if flag for details if T, print details about any resut
    write("RRegr package | eNanoMapper",                  file=outFile,append=T)
    write.table(paste("Regression method: ", RegrMethod), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("Split no.: ", iSplit), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("CV type: ", sCV),      file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Training Set Summary: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(summary(my.datf.train),   file=outFile,append=T,sep=",",col.names=T,quote=F)
    write.table("Test Set Summary: ",  file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(summary(my.datf.test), file=outFile,append=T,sep=",",col.names=T,quote=F)   
    
    write.table("Fitting Summary: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(data.frame(summary(glm.fit)$coefficients), file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Predictors: ",      file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(glm.fit), file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Trainig Results: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(lm.train.res), file=outFile,append=T,sep=",",col.names=T,quote=F)
    write.table("Test Results: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(lm.test.res), file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Full Statistics: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(my.stats, file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    # Variable Importance (max top 20)
    FeatImp <- varImp(glm.fit, scale = F)
    components = length(FeatImp)  # default plot all feature importance
    if (length(FeatImp)>20){     # if the number of features is greater than 20, use only 20
      components = 20
    }
    # Append feature importance to output details
    AppendList2CSv(FeatImp,outFile)
    
    fitModel <- glm.fit$finalModel
    
    # =============================================================================
    # Assessment of Applicability Domain (plot leverage)
    # =============================================================================
    
    # Residuals
    resids <- residuals(fitModel) # residuals
    write.table("Residuals of the fitted model: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(data.frame(resids), file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write residuals
    
    # Leverage / Hat values
    hat.fit <- hatvalues(fitModel)          # hat values
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
    
    # Cook's distance
    cook.dists<- cooks.distance(fitModel)
    cutoff.Cook <- 4/((nrow(my.datf.train)-length(fitModel$coefficients)-2)) # Cook's distance cutoff
    
    write.table("Cook's distances output: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("Cook's distance cutoff: ", cutoff.Cook), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Cook's distances: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(data.frame(cook.dists), file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write residuals
    
    # Influence
    infl <- influence(fitModel)#produces several statistics of the kind
    
    write.table("Point influence output: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Influences: ",file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(data.frame(infl), file=outFile,append=T,sep=",",col.names=T,row.names=T, quote=F) # write residuals
    
    # PDF with 12 plots
    # --------------------------------------------------------------
    pdf(file=paste(outFile,".",sCV,".","split",iSplit,".pdf",sep=""))
    # par(mfrow = c(3, 4)) # all plots into one page!
    
    plot(my.datf.train[,1],pred.tr,xlab="Yobs", ylab="Ypred", type="b", main="Train Yobs-Ypred") # plot 1
    plot(my.datf.test[,1], pred.ts,xlab="Yobs", ylab="Ypred", type="b", main="Test Yobs-Ypred")  # plot 2
    if(length(is.na(c(FeatImp$importance$Overall)))<=(length(c(FeatImp$importance$Overall))-3)){  
      dotchart(as.matrix(FeatImp$importance),main="Feature Importance")}                          # plot 3
    
    # Fitted vs Residuals - plot 4
    plot(fitted(fitModel),residuals(fitModel),
         main="Fitted vs. Residuals for Fitted Model",
         xlab="Fitted", ylab="Residuals")
    abline(h = 0, lty = 2)
    
    # Leverage plots - plot 5
    plot(hat.fit, type = "h",
         main="Leverage for Fitted Model",
         xlab="Index", ylab="Hat")
    abline(h = thresh.lever, lty = 2, col="red") # leverage thresh
    
    # Cook's distance - plot 6
    if(length(is.na(cook.dists))<=(length(cook.dists)-3)){  
      plot(cook.dists,
           main="Cook's Distance for Fitted Model",
           xlab="Index", ylab="Cook Distance")
      
      for (p in 1:6) {
        plot(fitModel, which=p, cook.levels=cutoff.Cook) # 6 standard fitting plots
      }
    }    
    # plot(FeatImp, top = components,main="Feature Importance") # ERROR !
    dev.off()
    # --------------------------------------------------------------
  }
  return(list(stat.values= my.stats, model=glm.fit))  # return a list with statistics and the full model
}
#----------------------------------------------------------------------------------------------------------------------
