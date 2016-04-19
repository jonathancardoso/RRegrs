# -----------------------------------------------------------------------
# svm regression function helper
# -----------------------------------------------------------------------
# jseoane
# use:
# svmFuncsGradW: RAKOTOMAMONJY gradient w
load(system.file("models", "model.svmRadialReg.RData", package = "RRegrs"))

svmFuncsW = caretFuncs    ## regular ranking using w
svmFuncsW$fit=function(x,y,first,last,...,tuneGrid){
  #cat(param$sigma,"\n")
  #library(kernlab)
  sigma = sigest(x)[2]  
  cs = tuneGrid$.C
  eps = tuneGrid$.epsilon
  tuneGrid = expand.grid(.C=cs,.sigma=sigma,.epsilon=eps)  
  train(x,y,...,tuneGrid=tuneGrid)
}
#--------------------------------------------------------------------------------
svmFuncsW$rank=function(object,x,y){
  alphas = alpha(object$finalModel)
  alpha.idxs = alphaindex(object$finalModel)
  y.sv = as.numeric(y[alpha.idxs])
  w = (y.sv * alphas) %*% xmatrix(object$finalModel)
  sig = ifelse(object$finalModel@fitted>y,yes=1,no=-1)
  
  avImp = t(w*w)
  out = data.frame(avImp)
  colnames(out) = "Overall"
  out = out[order(out$Overall, decreasing = TRUE), , drop = FALSE]
  out$var <- rownames(out)
  out
}
#--------------------------------------------------------------------------------
svmFuncsW$pred= function(object, x)
{
  tmp = predict(object, newdata=x)
  if(object$modelType == "Classification" &
     !is.null(object$modelInfo$prob))
  {
    out1 =cbind(data.frame(pred = tmp),
                as.data.frame(predict(object$finalModel, newdata=x, type = "prob")))
  } else out1 <- tmp
  out1
}
#--------------------------------------------------------------------------------

# Based on the gradient of svm coefs
svmFuncsGradW = svmFuncsW
svmFuncsGradW$rank=function(object,x,y){ # RAKOTOMAMONJY gradient w  
  alphas = alpha(object$finalModel)#[[1]]
  alpha.idxs = alphaindex(object$finalModel)#[[1]]
  y.sv = y[alpha.idxs]  
  krnFun = kernelf(object$finalModel)
  kernel = kernelMatrix(krnFun,x)
  sigma = krnFun@kpar$sigma
  xmat = xmatrix(object$finalModel)[[1]]
  kerSV = kernel[alpha.idxs,alpha.idxs]  
  nSV = length(alpha.idxs)
  nfeat = dim(x)[2]
  avImp = numeric(nfeat)
  names(avImp) = colnames(x)
  
  for(i in 1:nfeat){
    deraux =  (  x[alpha.idxs,i] %*% t(as.matrix(rep(1,nSV))) ) -  (as.matrix(rep(1,nSV)) %*% t(x[alpha.idxs,i])     )    
    kernelDeriv1 = -(deraux * kerSV) / (sigma^2)
    kernelDeriv2 =  (deraux * kerSV) / (sigma^2)
    gradMarg1= -t(y.sv*alphas) %*% kernelDeriv1 %*%  (y.sv*alphas)
    gradMarg2= -t(y.sv*alphas) %*% kernelDeriv2 %*%  (y.sv*alphas)
    avImp[i] = gradMarg1^2 + gradMarg2^2
  }
  
  out = data.frame(avImp)
  colnames(out) = "Overall"
  out = out[order(out$Overall, decreasing = TRUE), , drop = FALSE]
  out$var <- rownames(out)
  out
}
#----------------------------------------------------------------------------------------------------------------------



SVRMreg <- function(my.datf.train,my.datf.test,sCV,iSplit=1,fDet=F,outFile="",cs=c(1,5,10,15,20)) {
  #====================================
  # 8.6 SVM Radial Regression (caret)
  #====================================
  
  #library(caret)
  #library(kernlab)
  
  cs = as.numeric(cs)
  
  net.c = my.datf.train[,1] # dependent variable is the first column in Training set
  RegrMethod <- "svmRadial" # type of regression
  
  # Define the CV conditions
  ctrl<- trainControl(method=sCV,number=10,repeats=10,
                      summaryFunction=defaultSummary)
  
  # Train the model using only training set
  set.seed(iSplit)
  sigma = sigest (as.matrix(my.datf.train[,-1]))[2]
  svmL.fit<- train(net.c~.,data=my.datf.train,
                   method='svmRadial',tuneLength=10,trControl=ctrl,
                   metric='RMSE',
                   tuneGrid=expand.grid(.sigma=sigma,.C= cs))
  
  #------------------------------
  # Training RESULTS
  #------------------------------
  RMSE.tr  <- svmL.fit$results[,which(names(svmL.fit$results)=='RMSE')]#2]
  R2.tr    <- svmL.fit$results[,which(names(svmL.fit$results)=='Rsquared')]#3]
  if (sCV == "repeatedcv"){ # if 10-fold CV
    RMSEsd.tr <- svmL.fit$results[,which(names(svmL.fit$results)=='RMSESD')]##4]
    R2sd.tr   <- svmL.fit$results[,which(names(svmL.fit$results)=='RsquaredSD')]#5]
  }
  if (sCV == "LOOCV"){ # if LOOCV
    RMSEsd.tr <- 0 # formulas will be added later!
    R2sd.tr   <- 0 # formulas will be added later!
  }
  
  #------------------------------------------------
  # RMSE & R^2, for train/test respectively
  #------------------------------------------------
  lm.train.res <- if(sCV != "none") { getTrainPerf(svmL.fit)}
  lm.test.res  <- postResample(predict(svmL.fit,my.datf.test),my.datf.test[,1])
  
  #------------------------------------------------
  # Adj R2, Pearson correlation
  #------------------------------------------------
  pred.tr     <- predict(svmL.fit,my.datf.train) # predicted Y
  pred.ts     <- predict(svmL.fit,my.datf.test)  # predicted Y
  noFeats.fit <- length(predictors(svmL.fit))    # no. of features from the fitted model
  Feats.fit   <- paste(predictors(svmL.fit),collapse="+") # string with the features included in the fitted model
  
  ds.full     <- rbind(my.datf.train,my.datf.test)
  pred.both   <- predict(svmL.fit,ds.full)       # predicted Y
  mae.tr      <- mae(my.datf.train[,1],pred.tr)
  mae.ts      <- mae(my.datf.test[,1],pred.ts)
  adjR2.tr    <- r2.adj.funct(my.datf.train[,1],pred.tr,noFeats.fit)
  adjR2.ts    <- r2.adj.funct(my.datf.test[,1],pred.ts,noFeats.fit)
  corP.ts     <- cor(my.datf.test[,1],pred.ts)
  
  adjR2.both  <- r2.adj.funct(ds.full[,1],pred.both,noFeats.fit)
  RMSE.both   <- rmse.funct(ds.full[,1],pred.both)
  r2.both     <- r2.funct(ds.full[,1],pred.both)
  
  # Generate the output list with statistics for each cross-validation type
  # ------------------------------------------------------------------------
  my.stats <- list("RegrMeth"     = RegrMethod,
                   "Split No"     = as.numeric(iSplit),     # from function param
                   "CVtype"       = sCV,                    # from function param
                   "NoModelFeats" = as.numeric(noFeats.fit),
                   "ModelFeats"   = Feats.fit,
                   "adjR2.tr"  = as.numeric(adjR2.tr),
                   
                   "RMSE.tr"   = as.numeric(min(RMSE.tr)),  # these 4 lines correspond to the min of RMSE.tr !!!
                   "MAE.tr"    = as.numeric(mae.tr),
                   "R2.tr"     = as.numeric(R2.tr[which.min(RMSE.tr)]),  
                   "RMSEsd.tr" = as.numeric(RMSEsd.tr[which.min(RMSE.tr)]),
                   "R2sd.tr"   = as.numeric(R2sd.tr[which.min(RMSE.tr)]),
                   
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
  if (fDet==T) {   # if flag for details if true, print details about any resut
    write("RRegr package | eNanoMapper", file=outFile,append=T)
    write.table(paste("Regression method: ", RegrMethod), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("Split no.: ", iSplit), file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(paste("CV type: ", sCV),      file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table("Training Set Summary: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(summary(my.datf.train),   file=outFile,append=T,sep=",",col.names=T,quote=F)
    write.table("Test Set Summary: ",  file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(summary(my.datf.test), file=outFile,append=T,sep=",",col.names=T,quote=F)   
    
    
    write.table("Predictors: ",       file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(svmL.fit), file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Trainig Results: ",      file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(lm.train.res), file=outFile,append=T,sep=",",col.names=T,quote=F)
    write.table("Test Results: ",        file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(predictors(lm.test.res), file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Full Statistics: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(my.stats, file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    # Variable Importance (max top 20)
    FeatImp <- varImp(svmL.fit, scale = F)
    components = length(FeatImp)  # default plot all feature importance
    if (length(FeatImp)>20){     # if the number of features is greater than 20, use only 20
      components = 20
    }
    # Append feature importance to output details
    AppendList2CSv(FeatImp,outFile)
    
    fitModel <- svmL.fit$finalModel
    
    # =============================================================================
    # Assessment of Applicability Domain (plot leverage)
    # =============================================================================
    
    # Residuals    
    resids = pred.tr-svmL.fit$trainingData$.outcome # residuals
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
      
      # Cook's distance ?
    }
    
    # Influence ?
    
    # PDF plots
    # --------------------------------------------------------------
    pdf(file=paste(outFile,".",sCV,".","split",iSplit,".pdf",sep=""))
    plot(my.datf.train[,1],pred.tr,xlab="Yobs", ylab="Ypred", type="b", main="Train Yobs-Ypred")
    plot(my.datf.test[,1], pred.ts,xlab="Yobs", ylab="Ypred", type="b", main="Test Yobs-Ypred")
    dotchart(as.matrix(FeatImp$importance),main="Feature Importance")
    
    # Fitted vs Residuals
    plot(fitted(fitModel),resids,
         main="Fitted vs. Residuals for Fitted Model",
         xlab="Fitted", ylab="Residuals")
    abline(h = 0, lty = 2)
    
    # Leverage plots
    if(det.Traind.pls!=0){
      plot(hat.fit, type = "h",
           main="Leverage for Fitted Model",
           xlab="Index", ylab="Hat")
      abline(h = thresh.lever, lty = 2, col="red") # leverage thresh
    }
    
    dev.off()
    # --------------------------------------------------------------
  }
  
  return(list(stat.values=my.stats, model=svmL.fit))  # return a list with statistics and the full model
}
#----------------------------------------------------------------------------------------------------------------------

SVMRFEreg <- function(my.datf.train,my.datf.test,sCV,iSplit=1,fDet=F,outFile="",cs=c(1,5,10,15,20),eps=c(0.01,0.1,0.3)) {
  #SVMRFEreg <- function(my.datf.train,my.datf.test,sCV,iSplit=1,fDet=F,outFile="",cs=c(1:10),eps=c(0.01,0.1,0.3),noCores=1) {
  
  #===========================================
  # SVM-RFE
  #===========================================
  
  #library(kernlab)
  
  net.c = my.datf.train[,1]   # make available the names of variables from training dataset
  RegrMethod <- "svmRFE" # type of regression
  
  # Define the CV conditions
  ctrl<- trainControl(method=sCV, number=3,repeats=1,#number=10,repeats=10,
                      summaryFunction=defaultSummary,verboseIter = F)
  
  rfeCtr = rfeControl(functions=svmFuncsGradW,method="cv",number=10,repeats=10, saveDetails = T, verbose=T,rerank = T,allowParallel=T)#number=10,repeats=10,
  sigma = sigest (as.matrix(my.datf.train[,-1]))[2]
  #cs = c(0.0001,0.1,1,5,15,50)  # EXTERNAL PARAMETERS!!!
  #cs = c(1,5,15,50)
  #eps=c(0.01,0.1,0.3)   # EXTERNAL PARAMETERS!!!
  sizes = 2^(1:sqrt(ncol(my.datf.train)-1))
  tuneVars = expand.grid(.C=cs, .sigma=sigma, .epsilon=eps)    
  
  # Train the model using only training set
  set.seed(iSplit)
  
  rfesvm.fit = rfe(as.matrix(my.datf.train[,-1]),net.c,sizes = sizes,rfeControl=rfeCtr,prob.model =F,method=svmRadialReg,tuneGrid = tuneVars,trControl=ctrl ,allowParallel=T)
  
  # warning in some of the parameters is a extreme value
  if(rfesvm.fit$fit$bestTune$C %in% cs[c(1,length(cs))])
    warning("Best fitted value of C=",rfesvm.fit$fit$bestTune$C," is a extreme value in your possible c values. You may want to reset your C paramenter options", call. = FALSE)
  if(rfesvm.fit$fit$bestTune$epsilon %in% eps[c(1,length(cs))])
    warning("Best fitted value of eps=",rfesvm.fit$fit$bestTune$epsilon," is a extreme value in your possible epsilon values. You may want to reset your eps paramenter options", call. = FALSE)
  
  #------------------------------
  # Training RESULTS
  #------------------------------
  RMSE.tr  <- rfesvm.fit$results[rfesvm.fit$results$Variables ==rfesvm.fit$bestSubset,2]
  R2.tr    <- rfesvm.fit$results[rfesvm.fit$results$Variables ==rfesvm.fit$bestSubset,3]
  if (sCV == "repeatedcv"){ # if 10-fold CV
    RMSEsd.tr <- rfesvm.fit$results[rfesvm.fit$results$Variables ==rfesvm.fit$bestSubset,4]
    R2sd.tr   <- rfesvm.fit$results[rfesvm.fit$results$Variables ==rfesvm.fit$bestSubset,5]
  }
  if (sCV == "LOOCV"){ # if LOOCV
    RMSEsd.tr <- 0 # formulas will be added later!  TODOOOOOOOOOOOOOO
    R2sd.tr   <- 0 # formulas will be added later!
  }
  
  #------------------------------------------------
  # RMSE & R^2, for train/test respectively
  #------------------------------------------------
  rfesvm.train.res <- getTrainPerf(rfesvm.fit$fit)
  rfesvm.test.res  <- postResample(predict(rfesvm.fit,my.datf.test),my.datf.test[,1])
  
  #------------------------------------------------
  # Adj R2, Pearson correlation
  #------------------------------------------------
  pred.tr     <- predict(rfesvm.fit,my.datf.train) # predicted Y for training
  pred.ts     <- predict(rfesvm.fit,my.datf.test)  # predicted Y for test
  noFeats.fit <- length(predictors(rfesvm.fit))    # no. of features from the fitted model
  Feats.fit   <- paste(predictors(rfesvm.fit),collapse="+") # string with the features included in the fitted model
  
  ds.full     <- rbind(my.datf.train,my.datf.test)
  pred.both   <- predict(rfesvm.fit,ds.full)       # predicted Y
  mae.tr      <- mae(my.datf.train[,1],pred.tr)
  mae.ts      <- mae(my.datf.test[,1],pred.ts)
  adjR2.tr    <- r2.adj.funct(my.datf.train[,1],pred.tr,noFeats.fit)
  adjR2.ts    <- r2.adj.funct(my.datf.test[,1],pred.ts,noFeats.fit)
  corP.ts     <- cor(my.datf.test[,1],pred.ts)
  
  adjR2.both  <- r2.adj.funct(ds.full[,1],pred.both,noFeats.fit)
  RMSE.both   <- rmse.funct(ds.full[,1],pred.both)
  r2.both     <- r2.funct(ds.full[,1],pred.both)
  
  # Generate the output list with statistics for each cross-validation type
  # -----------------------------------------------------------------------------
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
    write.table(predictors(rfesvm.fit), file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Trainig Results: ",     file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(rfesvm.train.res,file=outFile,append=T,sep=",",col.names=T,quote=F)
    write.table("Test Results: ",        file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(rfesvm.test.res, file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    write.table("Full Statistics: ", file=outFile,append=T,sep=",",col.names=F,row.names=F,quote=F)
    write.table(my.stats,            file=outFile,append=T,sep=",",col.names=T,quote=F)
    
    # Variable Importance (max top 20)
    
    FeatImp <- svmFuncsGradW$rank(rfesvm.fit$fit,as.matrix(ds.full[,-1]),ds.full[,1])
    FeatImp = FeatImp[order(FeatImp[,1],decreasing=T),]
    components = length(FeatImp)  # default plot all feature importance
    if (length(FeatImp)>20){     # if the number of features is greater than 20, use only 20
      components = 20
    }
    # Append feature importance to output details
    AppendList2CSv(FeatImp,outFile)
    
    fitModel <- rfesvm.fit$fit
    
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
    fi = as.matrix(FeatImp[,1])
    rownames(fi)=FeatImp[,2]
    dotchart(fi,main="Feature Importance") # plot 3
    
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
  
  return(list(stat.values=my.stats, model=rfesvm.fit))  # return a list with statistics and the full model
}
#----------------------------------------------------------------------------------------------------------------------

