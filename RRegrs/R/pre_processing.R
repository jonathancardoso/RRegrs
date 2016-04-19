# ************************************
# RRegrs Specific functions
# ************************************

RemNear0VarCols <- function(ds,fDet=FALSE,outFile="ds3.No0Var.csv") {
  #================================================
  # Removal of near zero variance columns (Step 3)
  #================================================
  # inputs:
  # - ds = dataset frame
  # - fDet = flag for detais (TRUE/FALSE)
  # - outFileName = new file name  (it could include the path)
  
  # output = ds.Rem0NearVar  (ds without columns with near zero variance)
  # if datails = TRUE, output the new ds as a file
  # ------------------------------------------
  
  # default parameters are no details, with a CSV file name
  #library(caret)
  ds.Rem0NearVar <- ds               # default output without any modification
  ds.var <- nearZeroVar(ds)          # get the near zero columns
  if (!length(ds.var) == FALSE) {    # remove the columns only if nearZeroVar identified; if no columns to remove, ds will be the same
    ds.Rem0NearVar <- ds[,-(ds.var)] # get only the columns without this problem
    if (fDet == TRUE) {              # write as details the corrected ds file
      write.csv(ds.Rem0NearVar, outFile,row.names=F, quote=F)
    }
  }
  return(as.data.frame(ds.Rem0NearVar)) # return the new data frame without near zero variance
}
#----------------------------------------------------------------------------------------------------------------------

ScalingDS <- function(ds,s=1,c=2,fDet=FALSE,outFileName="ds4.scaled.csv") {
  #===========================
  # Scaling dataset (Step 4)
  #===========================
  # s = { 1,2,3 } - type of scaling: 1 = normalization, 2 = standardization, 3 = other
  # c = the number of column into the dataset to start scaling
  # - if c = 1: included the dependent variable
  # - if c = 2: only the features will be scaled
  # fDet = if details need to be printed    (TRUE/FALSE)
  # outFileName = new file name  (it could include the path)
  # Default scaling = NORMALIZATION !
  
  # DEFAULT scaled dataset = original
  # if other s diffent of 1,2,3 is used => no scaling!
  DataSet.scaled <- ds
  
  # if NORMALIZATION
  if (s==1) {
    # Scale all the features (from column c; column 1 is the predictor output)
    if(c==2){
      maxs <- apply(ds[c:ncol(ds)], 2, max)
      mins <- apply(ds[c:ncol(ds)], 2, min)
      ds.norm.scale<-scale(ds[c:ncol(ds)], center = mins, scale = maxs - mins)
      DataSet.scaled<-cbind(ds[,1],ds.norm.scale)
    }else{
      maxs <- apply(ds, 2, max)
      mins <- apply(ds, 2, min)
      DataSet.scaled<-scale(ds, center = mins, scale = maxs - mins)
    }
    
  }
  # if STADARDIZATION
  if (s==2) {
    # Scale all the features (from column c; column 1 is the predictor output)
    if(c==2){
      DataSet.scaled <- scale(ds[c:ncol(ds)],center=TRUE,scale=TRUE)
      DataSet.scaled<-cbind(ds[,1],DataSet.scaled)
    }else{
      DataSet.scaled<-scale(ds,center=TRUE,scale=TRUE)
    }
  }
  
  # if other scaling
  if (s==3) {
    # Scale all the features (from feature 2 bacause feature 1 is the predictor output)
    # TO ADD THE CODE ! 
  }
  
  # if DETAILS
  if (fDet ==TRUE) {
    # write the result into a separated file
    write.csv(DataSet.scaled, outFileName,row.names=F, quote=F)  
  }
  return (as.data.frame(DataSet.scaled)) # return the scaled data frame
}
#----------------------------------------------------------------------------------------------------------------------

RemCorrs <- function(ds,fDet,cutoff,outFile) {
  # ========================================
  # Remove the correlated columns (Step 5)
  # ========================================
  # ds = dataset frame
  # fDet = flag fro details (TRUE/FALSE)
  # cutoff = correlation cut off (ex: 0.9)
  # outFileName = new file name  (it could include the path)
  
  # Generates 5 file results:
  # - returns a dataset without the correlated columns  (1 file)
  # - generate initial correlation matrix 
  #   and the one after removing the correlated features (2 files)
  # - plots for the before and after correlation removal (2 files)
  # ------------------------------------------------------------------------
  
  # another version of this function should be implemented using
  # pairwise test between i and j descriptors- if(r2>=0.9){remove the j descriptor}
  # using findCorelations() from caret
  
  #library(corrplot) #corrplot: the library to compute correlation matrix.
  #library(caret)
  
  DataSet <- ds               # input dataset
  DataSetFiltered.scale <- ds # default results without any modification
  
  # calculate the correlation matrix for the entire file!
  # !!! NEED TO BE CORRECTED to avoid dependent variable (first column) but to report it!
  corrMat <- cor(DataSet)                                             # get corralation matrix
  
  if (fDet==TRUE) {
    CorrMatFile <- paste(outFile,".corrMAT.csv",sep='')
    # write correlation matrix as output file
    write.csv(corrMat, CorrMatFile, row.names=F, quote=F)
    
    # Plot the matrix, clustering features by correlation index
    # corrplot(corrMat, order = "hclust")
    
    # plot the correlatio plot before correlation removal
    CorrPlotFile <-  paste(outFile,".corrs.png",sep='')
    png(height=1200, width=1200, pointsize=25, filename=CorrPlotFile)
    col1 <-rainbow(100, s = 1, v = 1, start = 0, end = 0.9, alpha = 1)
    corrplot(corrMat,tl.cex=3,title="Initial feature correlation matrix",
             method="circle",is.corr=FALSE,#type="full",
             cl.lim=c(-1,1),cl.cex=2,addgrid.col="red",
             addshade="positive",col=col1,
             addCoef.col = rgb(0,0,0, alpha = 0.6), mar=c(0,0,1,0), diag= FALSE) 
    dev.off()
  }
  
  highlyCor <- findCorrelation(corrMat, cutoff) # find corralated columns
  
  # if no correlation found, return the original dataset
  if (length(highlyCor) == 0){
    return (ds)
  }
  
  # Apply correlation filter with the cutoff only if exists!
  # by removing all the variable correlated with more than cutoff
  DataSetFiltered.scale <- DataSet[,-highlyCor]
  
  if (fDet==TRUE) {
    corrMat <- cor(DataSetFiltered.scale)
    # plot again the rest of correlations after removing the correlated columns
    #corrplot(corrMat, order = "hclust")
    
    # plot the correlation plot AFTER correlation removal
    #CorrPlotFile2 =  paste(outFile,".afterRemCorr.png",sep='')
    #png(height=1200, width=1200, pointsize=25, file=CorrPlotFile2)
    #col1 <-rainbow(100, s = 1, v = 1, start = 0, end = 0.9, alpha = 1)
    #corrplot(corrMat,tl.cex=3,title="Correlation matrix after removing correlated features",
    #         method="circle",is.corr=FALSE,type="full",
    #         cl.lim=c(-1,1),cl.cex=2,addgrid.col="red",
    #         addshade="positive",col=col1,
    #         addCoef.col = rgb(0,0,0, alpha = 0.6), mar=c(0,0,1,0), diag= FALSE) 
    #dev.off()
    # correlation matrix for the rest of the columns after removal
    #CorrMatFile2 <- paste(outFile,".corrMAT4Selected.csv",sep='')
    # write correlation matrix as output file
    #write.csv(corrMat, CorrMatFile2, row.names=F, quote=F)
    # write the new dataset without the correlated features
    write.csv(DataSetFiltered.scale, outFile, row.names=F, quote=F)
  }
  
  return(as.data.frame(DataSetFiltered.scale))
}
#----------------------------------------------------------------------------------------------------------------------