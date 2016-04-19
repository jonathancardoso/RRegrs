# ======================================================================
# RRegrs - R Regressions
# ======================================================================
# Get the best regression models for one dataset using R caret methods
# eNanoMapper.net
# -------------------------------------------------------------------------------------------------------------
# AUTHORS: 
# -------------------------------------------------------------------------------------------------------------
# Georgia Tsiliki: ChemEng - NTUA, Greece, g_tsiliki@hotmail.com
# Cristian R. Munteanu: RNASA-IMEDIR, University of A Coruna, Spain, muntisa@gmail.com
# Jose A. Seoane: Stanford Cancer Institute, USA, seoane@stanford.edu
# Carlos Fernandez-Lozano: RNASA-IMEDIR, University of A Coruna, Spain, carlos.fernandez@udc.es
# Haralambos Sarimveis: ChemEng - NTUA, Greece, hsarimv@central.ntua.gr
# Egon Willighagen: BiGCaT - Maastricht University, Netherlands, egon.willighagen@gmail.com
# -------------------------------------------------------------------------------------------------------------


###############################################################################################
# RRegrs MAIN FUNCTION 
###############################################################################################

RRegrs <- function(DataFileName="ds.House.csv",DataFileSep=",",PathDataSet="DataResults",noCores=1,
  ResAvgs="RRegsResAvgs.csv",ResBySplits="RRegrsResAllSplits.csv",ResBest="RRegrsResBest.csv",
  fDet="T",fFilters="F",fScaling="T",fRemNear0Var="T",fRemCorr="T",
  fLM="T",fGLM="T",fPLS="T",fLASSO="T",fSVRM="T",fNN="T",fRF="T",fRFRFE="T",fSVMRFE="T",fENET="T",
  RFE_SVM_C="1;5;10;15;20",RFE_SVM_epsilon="0.01;0.1;0.3",
  cutoff=0.9,iScaling=1,iScalCol=1,trainFrac=0.75,iSplitTimes=10,noYrand=100,
  CVtypes="repeatedcv;LOOCV",NoNAValFile="ds.NoNA.csv",
  No0NearVarFile="ds.No0Var.csv",ScaledFile="ds.scaled.csv",NoCorrFile="ds.scaled.NoCorrs.csv",
  lmFile="LM.details.csv",glmFile="GLM.details.csv",plsFile="PLS.details.csv",
  lassoFile="Lasso.details.csv",svrmFile="SVMRadial.details.csv",
  nnFile="NN.details.csv",rfFile="RF.details.csv",rfrfeFile="RFRFE.details.csv",svmrfeFile="SVMRFE.details.csv",
  enetFile="ENET.details.csv",fR2rule="T") { # input = file with all parameters

  methodCount = 0;
  if (fLM=="T") methodCount = methodCount + 1;
  if (fGLM=="T") methodCount = methodCount + 1;
  if (fPLS=="T") methodCount = methodCount + 1;
  if (fLASSO=="T") methodCount = methodCount + 1;
  if (fSVRM=="T") methodCount = methodCount + 1;
  if (fNN=="T") methodCount = methodCount + 1;
  if (fRF=="T") methodCount = methodCount + 1;
  if (fRFRFE=="T") methodCount = methodCount + 1;
  if (fSVMRFE=="T") methodCount = methodCount + 1;
  if (fENET=="T") methodCount = methodCount + 1;

  if (methodCount < 2) stop("You must select at least two modelling methods to compare.");
  
  # fRBFdda="T", rbfDDAFile="RBF_DDA.details.csv",negThrStep=0.5
  
  # Minimal use:
  # RRegrs()                              # all default params
  # RRegrs(DataFileName="MyDataSet.csv")
  # RRegrs(DataFileName="MyDataSet.csv",PathDataSet="MyResultsFolder")
  # Default: all methods, no feature selection
  ptmTot <- proc.time() # total time

  # ----------------------------------
  # Parallel support
  # ----------------------------------
  if (noCores==0 | noCores>1){ # all available CPU cores or specific no of cores (if noCores = 1, no parallel support!)
    #noCoresSys=as.numeric(Sys.getenv('NUMBER_OF_PROCESSORS')) # automatically detected no. of CPU cores
    #library(parallel)
    noCoresSys=detectCores()
    
    if (noCores==0 | noCores>noCoresSys){ # all available CPU cores or the specific cores is greater than the available ones
      noCores=noCoresSys # use the available no of cores
    }

    # parallel for Linux or Mac:
    # ------------------------------------------
    if ( Sys.info() [['sysname']] == "Linux" | Sys.info() [['sysname']] == "Darwin" ){
      #library(doMC)
    }
    # ------------------------------------------
    # parallel for windows:
    # ------------------------------------------
    if (Sys.info()[['sysname']]=="Windows"){
      #library(doSNOW)
      #library(foreach)
    }  
  }     
  
  #==========================================================================================
  # (1) Load dataset and parameters
  #==========================================================================================
  
  # (1.1) PARAMETERS
  #------------------------------------------
  # Write parameter file
  #------------------------------------------
  dir.create(PathDataSet, showWarnings = FALSE)
  ParamFile <- file.path(PathDataSet, "Parameters.csv") # file to output the parameters

  # define a data frame with all parameters of the current calculation
  Params.df = data.frame(RRegrs.Parameters="DataFileName",Parameter.Value=as.character(DataFileName),Description="Input dataset file (Step 1)") # data frame with used parameters
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="PathDataSet",Parameter.Value=as.character(PathDataSet),Description="Working folder for all input and output files"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="noCores",Parameter.Value=as.character(noCores),Description="No of CPU cores (0=all available; 1=no parallel; >1 = specific no. of cores)"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="ResAvgs",Parameter.Value=as.character(ResAvgs),Description="Output file averaged statistics (by splits) for each regression method"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="ResBySplits",Parameter.Value=as.character(ResBySplits),Description="Output file statistics for each splitting and each regression method"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="ResBest",Parameter.Value=as.character(ResBest),Description="Output file statistics for the best model"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="fDet",Parameter.Value=as.character(fDet),Description="If calculate and print details for all the functions"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="fFilters",Parameter.Value=as.character(fFilters),Description="If run Filters (Step 2)"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="fScaling",Parameter.Value=as.character(fScaling),Description="If Scaling dataset (Step 3)"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="fRemNear0Var",Parameter.Value=as.character(fRemNear0Var),Description="If run Removal of near zero variance columns (Step 4)"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="fRemCorr",Parameter.Value=as.character(fRemCorr),Description="If run Removal of correlated columns (Step 5)"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="fLM",Parameter.Value=as.character(fLM),Description="If run LM (Step 8.1)"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="fGLM",Parameter.Value=as.character(fGLM),Description="If run GLM (Step 8.2)"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="fPLS",Parameter.Value=as.character(fPLS),Description="If run PLS (Step 8.3)"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="fLASSO",Parameter.Value=as.character(fLASSO),Description="If run LASSO (Step 8.4)"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="fENET",Parameter.Value=as.character(fENET),Description="If run ENET (Step 8.5)"))
  # Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="fRBFdda",Parameter.Value=as.character(fRBFdda),Description="If run RBF DDA (Step 8.6)"))
  # Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="negThrStep",Parameter.Value=as.character(negThrStep),Description="Negative Threshold step parameter for RBF DDA (Step 8.6)"))

  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="fSVRM",Parameter.Value=as.character(fSVRM),Description="If run svmRadial.RMSE (Step 8.7)"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="fNN",Parameter.Value=as.character(fNN),Description="If run Neural Networks (Step 8.8)"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="fRF",Parameter.Value=as.character(fRF),Description="If run Random Forest (Step 8.9)"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="fRFRFE",Parameter.Value=as.character(fRFRFE),Description="If run Random Forest RFE (Step 8.10)"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="fSVMRFE",Parameter.Value=as.character(fSVMRFE),Description="If run Random Forest (Step 8.11)"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="RFE_SVM_C",Parameter.Value=as.character(RFE_SVM_C),Description="Values of C for SVM RFE"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="RFE_SVM_epsilon",Parameter.Value=as.character(RFE_SVM_epsilon),Description="Values of epsilon for SVM RFE"))
  
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="cutoff",Parameter.Value=as.character(cutoff),Description="Cut-off for correlated features (default = 0.9)"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="iScaling",Parameter.Value=as.character(iScaling),Description="Type of scaling: 1 = normalization; 2 = standardization; 3 = other; any other: no scaling"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="iScalCol",Parameter.Value=as.character(iScalCol),Description="Scaling columns: 1 = including dependent variable; 2: only all the features"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="trainFrac",Parameter.Value=as.character(trainFrac),Description="Fraction of training set from the entire dataset; the rest of dataset is the test set"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="iSplitTimes",Parameter.Value=as.character(iSplitTimes),Description="Number of splitting the dataset into train and test (default  = 10)"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="noYrand",Parameter.Value=as.character(noYrand),Description="Number of Y-Randomization (default = 100)"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="CVtypes",Parameter.Value=as.character(CVtypes),Description="Cross-validation types: 10-CV (repeatedcv) and LOOCV"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="NoNAValFile",Parameter.Value=as.character(NoNAValFile),Description="Dataset without NA values (if fDet is True)"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="No0NearVarFile",Parameter.Value=as.character(No0NearVarFile),Description="Dataset without zero near features from Step 3 (if fDet is True)"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="ScaledFile",Parameter.Value=as.character(ScaledFile),Description="Scaled dataset file from Step 4 (if fDet is True)"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="NoCorrFile",Parameter.Value=as.character(NoCorrFile),Description="Dataset after correction removal in Step 5 (if fDet is True)"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="lmFile",Parameter.Value=as.character(lmFile),Description="LM output file with details"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="glmFile",Parameter.Value=as.character(glmFile),Description="GLM output file with details"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="plsFile",Parameter.Value=as.character(plsFile),Description="PLS output file with details"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="lassoFile",Parameter.Value=as.character(lassoFile),Description="Lasso output file with details"))
  #Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="rbfDDAFile",Parameter.Value=as.character(rbfDDAFile),Description="RBF DDA output file with details"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="svrmFile",Parameter.Value=as.character(svrmFile),Description="SVM Radial output file with details"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="nnFile",Parameter.Value=as.character(nnFile),Description="NN output file with details"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="rfFile",Parameter.Value=as.character(rfFile),Description="RF output"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="rfrfeFile",Parameter.Value=as.character(rfrfeFile),Description="RF-RFE output"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="svmrfeFile",Parameter.Value=as.character(svmrfeFile),Description="SVM-RFE output"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="enetFile",Parameter.Value=as.character(enetFile),Description="ENET output"))
  Params.df = rbind(Params.df,data.frame(RRegrs.Parameters="fR2rule",Parameter.Value=as.character(fR2rule),Description="Best model rule: R2 (default = T) or adjR2 (F)"))

  write.csv(Params.df,file=ParamFile,row.names=F,quote=F) # write parameters to a CSV in the working folder

  # Get calculation parameters 
  fDet         = as.logical(fDet)         # flag to calculate and print details for all the functions
  fFilters     = as.logical(fFilters)     # flag to apply filters                          (2)
  fScaling     = as.logical(fScaling)     # flag for dataset Scaling                       (3)
  fRemNear0Var = as.logical(fRemNear0Var) # flag for Removal of near zero variance columns (4)
  fRemCorr     = as.logical(fRemCorr)     # flag for Removal of correlated columns         (5)
  
  cutoff       = as.numeric(as.character(cutoff))  # cut off for correlated features
  fLM          = as.logical(fLM)     # flag to run LM            (8.1)
  fGLM         = as.logical(fGLM)    # flag to run GLM           (8.2)
  fPLS         = as.logical(fPLS)    # flag to run PLS           (8.3)
  fLASSO       = as.logical(fLASSO)  # flag to run LASSO         (8.4)
  fenet        = as.logical(fENET)   # flag to run ElasticNet    (8.5)
  #fRBFdda      = as.logical(fRBFdda) # flat to run RBF DDA      (8.6)
  fSVRM        = as.logical(fSVRM)   # flat to run svmRadial     (8.7)
  fNN          = as.logical(fNN)     # flat to run NN            (8.8)
  fRF          = as.logical(fRF)     # flag to run RandomForest  (8.9)
  fRFRFE       = as.logical(fRFRFE)  # flag to run RF-RFE        (8.10)
  fSVMRFE      = as.logical(fSVMRFE) # flag to run SVM RFE       (8.11)
  
  rfe_SVM_param_c   = strsplit(as.character(RFE_SVM_C),";")[[1]] # values of C for SVM RFE
  rfe_SVM_param_eps = strsplit(as.character(RFE_SVM_epsilon),";")[[1]] # values of epsilon for SVM RFE
  # negThrStep

  fR2rule = as.logical(fR2rule) # flag to decide order rule for the best model (True for R2 and False for adjR2, default = True)
  # ----------------------------------------------------------------------------------------
  trainFrac   = as.numeric(as.character(trainFrac))   # the fraction of training set from the entire dataset; trainFrac = the rest of dataset, the test set
 
  CVtypes = strsplit(as.character(CVtypes),";")[[1]] # types of cross-validation methods
  CVtypes2 = c("repeatedcv") # for complex methods we run only 10-fold CV even the user is using other parameters!
  
  # Generate path + file name = original dataset
  if (file.exists(DataFileName)) { # is it a full path already?
    inFile <- DataFileName
  } else {
    inFile <- file.path(PathDataSet, DataFileName)
  }
  
  sDescription=paste("=======================================================================================================",
    "RRegrs - R Regression Models",
    "Get the best regression models for one dataset using R caret methods", "eNanoMapper.net","AUTHORS:",
    "Georgia Tsiliki: ChemEng - NTUA, Greece, g_tsiliki@hotmail.com",
    "Cristian R. Munteanu: RNASA-IMEDIR, University of A Coruna, Spain, muntisa@gmail.com",
    "Jose A. Seoane: Stanford Cancer Institute, USA, seoane@stanford.edu",
    "Carlos Fernandez-Lozano: RNASA-IMEDIR, University of A Coruna, Spain, carlos.fernandez@udc.es",
    "Haralambos Sarimveis: ChemEng - NTUA, Greece, hsarimv@central.ntua.gr",
    "Egon Willighagen: BiGCaT - Maastricht University, The Netherlands, egon.willighagen@gmail.com",
    "=======================================================================================================",sep="\n")
  cat(sDescription) # print package header information

  # -----------------------------------
  # (1.2) Load the ORIGINAL DATASET
  # -----------------------------------
  cat("\n-> Loading original dataset ...\n")  # it can contain errors, correlations, near zero variance columns
  cat("      ---> ",inFile,"\n")
  ds.dat0 <- read.csv(inFile,header=T,sep=DataFileSep)          # original dataset frame
  
  # resolving the text to number errors for future calculations
  ds.indx<- colnames(ds.dat0)[2:dim(ds.dat0)[2]]                    # FEATURE names (no dependent variable)
  ds.dat1<- ds.dat0[1:dim(ds.dat0)[1],2:dim(ds.dat0)[2]]            # dataset as columns
  ds.dat1<- apply(ds.dat1,1,function(x)as.numeric(as.character(x))) # dataset as row vectors to be used with caret!!!
  
  # dependent variable
  net.c<- ds.dat0[,1]
  net.c<- as.numeric(as.character(net.c)) # values
  # full ds frame with training and test
  ds<- as.data.frame(cbind(net.c,t(ds.dat1)))
  
  #========================================================
  # (2) FILTERS
  #     (it will be implemented in the future versions)
  #========================================================
  # 2.1 Outlier removal
  # 2.2 Custom filter (percentage threshold)
  # 2.3 Processing of missing values - use of preProcess();
  #     caret employs knnImpute algorithm to impute values from a neighborhood of k
  # TO BE IMPLEMENTED

  # -----------------------------------------------------------------------
  # (2) Remove NA values
  # -----------------------------------------------------------------------
  if(length(which(is.na(ds)==TRUE))!=0){
    cat("-> Removal of NA values ...\n")
    outFile <- file.path(PathDataSet,NoNAValFile) # the same folder as input  
    
    # get the ds without NA values-- currently use the default function (mean) 
    ds <- impute.funct(ds)
    if (fDet == TRUE){   # write as details the corrected ds file
     write.csv(ds, outFile,row.names=F, quote=F)}
  }
 
  # -----------------------------------------------------------------------
  # (3) Remove near zero variance columns
  # -----------------------------------------------------------------------
  if (fRemNear0Var==T) {
    cat("-> Removal of near zero variance columns ...\n")
    outFile <- file.path(PathDataSet,No0NearVarFile) # the same folder as input  
    
    # get the ds without near zero cols 
    ds <- cbind("net.c" = ds[,1],RemNear0VarCols(ds[,2:dim(ds)[2]],fDet,outFile))
    # use df without Y (predicted values), reconstruct the ds
    # inputs: ds, flag for details, output file
  }
  
  # -----------------------------------------------------------------------
  # (4) Scaling dataset: normalization (default), standardization, other
  # -----------------------------------------------------------------------
  if (fScaling==T) {
    cat("-> Scaling original dataset ...\n")
    outFile <- file.path(PathDataSet,ScaledFile)       # the same folder as input
    
    # run fuction for scaling input dataset file
    ds <- ScalingDS(ds,iScaling,iScalCol,fDet,outFile)
    # use df without Y (predicted values), reconstruct the ds
    # inputs: ds, type of scaling, flag for details, starting column, output file
  }
  
  # -----------------------------------------------------------------------
  # (5) Remove correlated features
  # -----------------------------------------------------------------------
  if (fRemCorr==T) {    
    cat("-> Removing correlated features ...\n") 
    outFile <- file.path(PathDataSet,NoCorrFile)    # the same folder as the input
    
    # run function to remove the correlations between the features
    ds <- cbind("net.c" = ds[,1],RemCorrs(ds[,2:dim(ds)[2]],fDet,cutoff,outFile))
  }
  
  # Check data has at least 5 columns for meaningful analysis
  if(dim(ds)[2] < 3 || dim(ds)[1] < 3){
  print(c(dim(ds)))
  stop(paste("Your corrected data set has dimensions:", paste(as.character(dim(ds)),collapse=', '),". Try repeating analysis without filtering options.",sep=''))
  }

  # print no of CPU cores used for calculation
  # noCoresSys=as.numeric(Sys.getenv('NUMBER_OF_PROCESSORS')) # automatically detected no. of CPU cores
  #library(parallel)
  noCoresSys=detectCores()
  if (noCores==0){ cat("       -> CPU Cores = ",noCoresSys,"(only complex methods)\n") }
  else{ cat("                  -> CPU Cores = ",noCores,   "(only complex methods)\n") }
  
  #=========================================================================================================
  # Steps 6 - 8 will be repeated 10 times for reporting each result and average (iSplitTimes = 10, default)
  #=========================================================================================================
  
  # Initialize the list with the statistics results; the same HEADER as the function output
  dfRes <- list("RegrMeth"     = NULL,
                "Split No"     = NULL,    
                "CVtype"       = NULL,      
                "NoModelFeats" = NULL,
                "ModelFeats"   = NULL,
                "adjR2.tr"  = NULL,
                "RMSE.tr"   = NULL,
                "R2.tr"     = NULL,
                "RMSEsd.tr" = NULL,
                "R2sd.tr"   = NULL,
                "adjR2.ts"= NULL,
                "RMSE.ts" = NULL,
                "R2.ts"   = NULL,
                "corP.ts" = NULL,
                "adjR2.both" = NULL,
                "RMSE.both"  = NULL,
                "R2.both"    = NULL)
  #-------------------------------------------------------------------------------------------------
  for (i in 1:iSplitTimes) {   # Step splitting number = i
    # -----------------------------------------------------------------------
    # (6) Dataset split: Training and Test sets
    # -----------------------------------------------------------------------
    cat("-> Splitting dataset in Training and Test sets ...\n")
    cat(paste("--> Split No.",i,"from",iSplitTimes,"\n"))

    # Initialize the list with the statistical models for all types of CV; per iSplitTimes
    dfMod <- sapply(CVtypes,function(x) NULL)
    for(cv in 1:length(CVtypes)){class(dfMod[[cv]])<- 'list'; names(dfMod)[[cv]]<- CVtypes[cv]}
    mod.ind<- rep(1,length(CVtypes)) # dummy variable to indicate the index of each new dfMod entry (per CVtype) 
    
    iSeed=i                 # to reapeat the ds splitting, different values of seed will be used
    dsList  <- DsSplit(ds,trainFrac,fDet,PathDataSet,iSeed) # return a list with 2 datasets = dsList$train, dsList$test
    
    # get train and test from the resulted list
    ds.train<- dsList$train
    ds.test <- dsList$test
    
    # -----------------------------------------------------------------------
    # (7) Feature selection
    # -----------------------------------------------------------------------
    # TO BE IMPLEMENTED
    
    # -----------------------------------------------------------------------
    # (8) REGRESSION METHODS
    # -----------------------------------------------------------------------
    
    # ----------------------------------
    # Parallel support - registed cores
    # ----------------------------------
    if (noCores!=1){
      # ------------------------------------------
      # parallel for Linux or Mac:
      # ------------------------------------------
      if (Sys.info()[['sysname']]=="Linux" | Sys.info()[['sysname']]=="Darwin"){
        registerDoMC(cores = noCores) # CPU cores
      }
      # ------------------------------------------
      # parallel for windows:
      # ------------------------------------------
      if (Sys.info()[['sysname']]=="Windows"){
       cl<-makeCluster(noCores,outfile="") 
        registerDoSNOW(cl)
      }   
    }
    
    # --------------------------------------------
    # 8.1. Basic LM : default
    # --------------------------------------------
    if (fLM==T) {   # if LM was selected, run the method
      outFile.LM <- file.path(PathDataSet,lmFile)   # the same folder as the input is used for the output
      cat("-> LM : Linear Multi-regression ...\n")
      
      # For each type of CV do all the statistics
      # -----------------------------------------------------
      for (cv in 1:length(CVtypes)) { # there is no CV but it will be implemented in the future!!!
        cat("    -->",CVtypes[cv],"\n")
        ptmLM <- proc.time()
        lm.model <- LMreg(ds.train,ds.test,CVtypes[cv],i,fDet,outFile.LM) # run GLM for each CV and regr method
        print(proc.time() - ptmLM) # print running time
        
        my.stats.LM <- lm.model$stat.values # stat values
        my.model.LM <- lm.model$model # model 
        
        #-------------------------------------------------------
        # Add output from GLM to the list of results
        #-------------------------------------------------------
        # List of results for each splitting, CV type & regression method
        dfRes = mapply(c, my.stats.LM, dfRes, SIMPLIFY=F)
        
        # List of models for each splitting, CV type & regression method
        names1 <- strsplit(deparse(quote(my.model.LM)),'my.model.')[[1]][2]
        dfMod[[cv]]$names1 <- my.model.LM  
        names(dfMod[[cv]])[mod.ind[cv]] <- names1[1]
        mod.ind[cv] <- mod.ind[cv] +1 # update mod.ind indicator variable
      } # end CV types
    } # end LM
    
    # -----------------------------------------------------------------------------------------
    # (8.2) GLM based on AIC regression - Generalized Linear Model with Stepwise Feature Selection
    # -----------------------------------------------------------------------------------------
    if (fGLM==T) {   # if GLM was selected, run the method
      outFile.GLM <- file.path(PathDataSet,glmFile)   # the same folder as the input is used for the output
      
      cat("-> GLM : Generalized Linear Model stepwise - based on AIC ...\n")
		  # For each type of CV do all the statistics
		  # -----------------------------------------------------
		  for (cv in 1:length(CVtypes)) {
		    cat("    -->",CVtypes[cv],"\n")
		    ptmGLM <- proc.time()
		    glm.model  <- GLMreg(ds.train,ds.test,CVtypes[cv],i,fDet,outFile.GLM) # run GLM for each CV and regr method
		    print(proc.time() - ptmGLM) # print running time
		    
		    my.stats.GLM <- glm.model$stat.values # stat values
		    my.model.GLM <- glm.model$model # model
		    #my.stats.split <- c(my.stats.dsInfo,my.stats.GLM) # merge the ds info with statistics results for each Cv & reg method
		    
		    #-------------------------------------------------------
		    # Add output from GLM to the list of results
		    #-------------------------------------------------------
		    # List of results for each splitting, CV type & regression method
		    dfRes = mapply(c, my.stats.GLM, dfRes, SIMPLIFY=F)
		    
		    # List of models for each splitting, CV type & regression method
		    names1 <- strsplit(deparse(quote(my.model.GLM)),'my.model.')[[1]][2]
		    dfMod[[cv]]$names1 <- my.model.GLM  
		    names(dfMod[[cv]])[mod.ind[cv]] <- names1[1]
		    mod.ind[cv] <- mod.ind[cv] +1 # update mod.ind indicator variable
		  } # end CV types  
    } # end GLM
    
    # --------------------------------------------
    # 8.3. PLS
    # --------------------------------------------
    if (fPLS==T) {   # if PLS was selected, run the method
		outFile.PLS <- file.path(PathDataSet,plsFile)   # the same folder as the input is used for the output
      
		cat("-> PLS : Partial Least Squares Regression ...\n")
		# For each type of CV do all the statistics
		# -----------------------------------------------------
		for (cv in 1:length(CVtypes)) {
			cat("    -->",CVtypes[cv],"\n")
			ptmPLS <- proc.time()
			pls.model <- PLSreg(ds.train,ds.test,CVtypes[cv],i,fDet,outFile.PLS) # run PLS for each CV and regr method
			print(proc.time() - ptmPLS) # print running time
        
			my.stats.PLS <- pls.model$stat.values
			my.model.PLS <- pls.model$model
			#-------------------------------------------------------
			# Add output from PLS to the list of results
			#-------------------------------------------------------
			# List of results for each splitting, CV type & regression method
			dfRes = mapply(c, my.stats.PLS, dfRes, SIMPLIFY=F)
        
			# List of models for each splitting, CV type & regression method
			names1 <- strsplit(deparse(quote(my.model.PLS)),'my.model.')[[1]][2]
			dfMod[[cv]]$names1 <- my.model.PLS  
			names(dfMod[[cv]])[mod.ind[cv]] <- names1[1]
			mod.ind[cv] <- mod.ind[cv] +1 # update mod.ind indicator variable
		} # end CV types
    } # end PLS
    
    #cat("-> PLS Wrapper Feature Selection ...\n")
    ## For each type of CV do all the statistics
    ## -----------------------------------------------------
    #for (cv in 1:length(CVtypes)) {
    #  cat("    -->",CVtypes[cv],"\n")
    #  ptmPLSw <- proc.time()
    #  pls.model <- PLSregWSel(ds.train,ds.test,CVtypes[cv],i,fDet,outFile.PLS) # run PLSw for each CV and regr method
    #  print(proc.time() - ptmPLSw) # print running time
    
    #  my.stats.PLS <- pls.model$stat.values
    #  my.model.PLS <- pls.model$model
    
    #  #-------------------------------------------------------
    #  # Add output from PLSw to the list of results
    #  #-------------------------------------------------------
    #  # List of results for each splitting, CV type & regression method
    #  dfRes = mapply(c, my.stats.PLS, dfRes, SIMPLIFY=F)
    #} # end CV types 
    
    # --------------------------------------------
    # 8.4. LASSO regression
    # --------------------------------------------
    if (fLASSO==T) {   # if LASSO was selected, run the method
		outFile.LASSO <- file.path(PathDataSet,lassoFile)   # the same folder as the input is used for the output
      
		cat("-> Lasso ...\n")
		# For each type of CV do all the statistics
		# -----------------------------------------------------
		for (cv in 1:length(CVtypes2)) {
			cat("    -->",CVtypes2[cv],"\n")
			ptmLASSO <- proc.time()
			lasso.model <- LASSOreg(ds.train,ds.test,CVtypes2[cv],i,fDet,outFile.LASSO) # run LASSO for each CV and regr method
			print(proc.time() - ptmLASSO) # print running time

			my.stats.LASSO <- lasso.model$stat.values
			my.model.LASSO <- lasso.model$model
			#-------------------------------------------------------
			# Add output from Lasso to the list of results
			#-------------------------------------------------------
			# List of results for each splitting, CV type & regression method
			dfRes = mapply(c, my.stats.LASSO, dfRes, SIMPLIFY=F)
			# List of models for each splitting, CV type & regression method
			names1 <- strsplit(deparse(quote(my.model.LASSO)),'my.model.')[[1]][2]
			dfMod[[cv]]$names1 <- my.model.LASSO  
			names(dfMod[[cv]])[mod.ind[cv]] <- names1[1]
			mod.ind[cv] <- mod.ind[cv] +1 # update mod.ind indicator variable
		} # end CV types
    } # end Lasso
    
    # --------------------------------------------
    # 8.5. Elastic Net regression
    # --------------------------------------------
    if (fenet==T) {   # if ENET was selected, run the method
		outFile.ENET <- file.path(PathDataSet,enetFile)   # the same folder as the input is used for the output
      
		cat("-> ENET : Elastic Nets ...\n")
		# For each type of CV do all the statistics
		# -----------------------------------------------------
		for (cv in 1:length(CVtypes)) {
			cat("    -->",CVtypes[cv],"\n")
			ptmENET <- proc.time()
			enet.model <- ENETreg(ds.train,ds.test,CVtypes[cv],i,fDet,outFile.ENET) # run elastic net for each CV and regr method
			print(proc.time() - ptmENET) # print running time
        
			my.stats.ENET <- enet.model$stat.values
			my.model.ENET <- enet.model$model 
        
			#-------------------------------------------------------
			# Add output from ENET to the list of results
			#-------------------------------------------------------
			# List of results for each splitting, CV type & regression method
			dfRes = mapply(c, my.stats.ENET, dfRes, SIMPLIFY=F)
			# List of models for each splitting, CV type & regression method
			names1 <- strsplit(deparse(quote(my.model.ENET)),'my.model.')[[1]][2]
			dfMod[[cv]]$names1 <- my.model.ENET  
			names(dfMod[[cv]])[mod.ind[cv]] <- names1[1]
			mod.ind[cv] <- mod.ind[cv] +1 # update mod.ind indicator variable
		} # end CV types  
    } # end enet

#     # ----------------------------------------------------------------
#     # 8.6. RBF network with the DDA algorithm regression (caret)
#     # ----------------------------------------------------------------
#     if (fRBFdda==T) {   # if RBF-DDA was selected, run the method
# 		outFile.rbfDDA <- file.path(PathDataSet,rbfDDAFile)   # the same folder as the input is used for the output
# 		  
# 		cat("-> RBF-DDA : Radial Basis Functions - Dynamic Decay Adjustment ...\n")
# 		# For each type of CV do all the statistics
# 		# -----------------------------------------------------
# 		for (cv in 1:length(CVtypes2)) {
# 		    cat("    -->",CVtypes2[cv],"\n")
# 		    ptmRBF_DDA <- proc.time()
# 		    rbfDDA.model <- RBF_DDAreg(ds.train,ds.test,CVtypes2[cv],negThrStep,i,fDet,outFile.rbfDDA) # run rbfDDA for each CV and regr method
# 		    print(proc.time() - ptmRBF_DDA) # print running time
# 		    
# 		    my.stats.rbfDDA <- rbfDDA.model$stat.values
# 		    my.model.rbfDDA <- rbfDDA.model$model 
# 		    #-------------------------------------------------------
# 		    # Add output from RBF-DDA to the list of results
# 		    #-------------------------------------------------------
# 		    # List of results for each splitting, CV type & regression method
# 		    dfRes = mapply(c, my.stats.rbfDDA, dfRes, SIMPLIFY=F)
# 		    
# 		    # List of models for each splitting, CV type & regression method
# 		    names1 <- strsplit(deparse(quote(my.model.rbfDDA)),'my.model.')[[1]][2]
# 		    dfMod[[cv]]$names1 <- my.model.rbfDDA  
# 		    names(dfMod[[cv]])[mod.ind[cv]] <- names1[1]
# 		    mod.ind[cv] <- mod.ind[cv] +1 # update mod.ind indicator variable
# 		} # end CV types
#     } # end rbfDDA
    
    # --------------------------------------------
    # 8.7. SVM radial regression
    # --------------------------------------------
    if (fSVRM==T) {   # if SVM Radial was selected, run the method
		outFile.SVRM <- file.path(PathDataSet,svrmFile)   # the same folder as the input is used for the output
      
		cat("-> SVM radial : Support vector machine using radial functions ...\n")
		# For each type of CV do all the statistics
		# -----------------------------------------------------
		for (cv in 1:length(CVtypes)) {
		    cat("    -->",CVtypes[cv],"\n")
		    ptmSVRM <- proc.time()
		    SVRM.model <- SVRMreg(ds.train,ds.test,CVtypes[cv],i,fDet,outFile.SVRM,rfe_SVM_param_c) # run SVRM Radial for each CV and regr method
		    print(proc.time() - ptmSVRM) # print running time
		    
		    my.stats.SVRM <- SVRM.model$stat.values
		    my.model.SVRM <- SVRM.model$model 
		    
		    #-------------------------------------------------------
		    # Add output from SVM Radial to the list of results
		    #-------------------------------------------------------
		    # List of results for each splitting, CV type & regression method
		    dfRes = mapply(c, my.stats.SVRM, dfRes, SIMPLIFY=F)
		    # List of models for each splitting, CV type & regression method
		    names1 <- strsplit(deparse(quote(my.model.SVRM)),'my.model.')[[1]][2]
		    dfMod[[cv]]$names1 <- my.model.SVRM  
		    names(dfMod[[cv]])[mod.ind[cv]] <- names1[1]
		    mod.ind[cv] <- mod.ind[cv] +1 # update mod.ind indicator variable
		} # end CV types
    } # end SVRM
    
    # --------------------------------------------
    # 8.8. Neural Networks Regression
    # --------------------------------------------
    if (fNN==T) {   # if NNet was selected, run the method
      outFile.NN <- file.path(PathDataSet,nnFile)   # the same folder as the input is used for the output
      
		 cat("-> NN : Neural Networks ...\n")
		 # For each type of CV do all the statistics
		 # -----------------------------------------------------
		 for (cv in 1:length(CVtypes)) {
			cat("    -->",CVtypes[cv],"\n")
		    ptmNN <- proc.time()
		    nn.model <- NNreg(ds.train,ds.test,CVtypes[cv],i,fDet,outFile.NN) # run NNet for each CV and regr method
		    print(proc.time() - ptmNN) # print running time
		    
		    my.stats.NN <- nn.model$stat.values
		    my.model.NN <- nn.model$model
		    
		    #-------------------------------------------------------
		    # Add output from NNet to the list of results
		    #-------------------------------------------------------
		    # List of results for each splitting, CV type & regression method
		    dfRes = mapply(c, my.stats.NN, dfRes, SIMPLIFY=F)
		    # List of models for each splitting, CV type & regression method
		    names1 <- strsplit(deparse(quote(my.model.NN)),'my.model.')[[1]][2]
		    dfMod[[cv]]$names1 <- my.model.NN  
		    names(dfMod[[cv]])[mod.ind[cv]] <- names1[1]
		    mod.ind[cv] <- mod.ind[cv] +1 # update mod.ind indicator variable
		} # end CV types  
    } # end NNet
    
    # --------------------------------------------
    # 8.9. Random Forest Regression (RF)
    # --------------------------------------------
    if (fRF==T) {   # if RF was selected, run the method
		outFile.RF <- file.path(PathDataSet,rfFile)   # the same folder as the input is used for the output
		cat("-> RF : Random Forest ...\n")
		# For each type of CV do all the statistics
		# -----------------------------------------------------
		for (cv in 1:length(CVtypes2)) {
			cat("    -->",CVtypes2[cv],"\n")
			ptmRF <- proc.time()
			rf.model <- RFreg(ds.train,ds.test,CVtypes[cv],i,fDet,outFile.RF) # run RF for each CV and regr method
			print(proc.time() - ptmRF) # print running time
		   
			my.stats.RF <- rf.model$stat.values
			my.model.RF <- rf.model$model
		  
			#-------------------------------------------------------
			# Add output from RF to the list of results
			#-------------------------------------------------------
			# List of results for each splitting, CV type & regression method
			dfRes = mapply(c, my.stats.RF, dfRes, SIMPLIFY=F)
			# List of models for each splitting, CV type & regression method
			names1 <- strsplit(deparse(quote(my.model.RF)),'my.model.')[[1]][2]
			dfMod[[cv]]$names1 <- my.model.RF  
			names(dfMod[[cv]])[mod.ind[cv]] <- names1[1]
			mod.ind[cv] <- mod.ind[cv] +1 # update mod.ind indicator variable
		} # end CV types
	}
      
	# --------------------------------------------
	# 8.10. RF-RFE: Random Forest Regression Recursive Feature Elimination
	# --------------------------------------------
	if (fRFRFE==T) {   # if RF-RFE was selected, run the method
		outFile.RFRFE <- file.path(PathDataSet,rfrfeFile)   # the same folder as the input is used for the output
      
		cat("-> RF-RFE: Random Forest-Recursive Feature Elimination ...\n")
		# For each type of CV do all the statistics
		# -----------------------------------------------------
		for (cv in 1:length(CVtypes2)) {
			cat("    -->",CVtypes2[cv],"\n")
			ptmRFRFE <- proc.time()
			rfrfe.model <- RFRFEreg(ds.train,ds.test,CVtypes[cv],i,fDet,outFile.RFRFE) # run RF for each CV and regr method
			print(proc.time() - ptmRFRFE) # print running time
        
			my.stats.RFRFE <- rfrfe.model$stat.values
			my.model.RFRFE <- rfrfe.model$model
			
			#-------------------------------------------------------
			# Add output from RF to the list of results
			#-------------------------------------------------------
			# List of results for each splitting, CV type & regression method
			dfRes = mapply(c, my.stats.RFRFE, dfRes, SIMPLIFY=F)
			# List of models for each splitting, CV type & regression method
			names1 <- strsplit(deparse(quote(my.model.RFRFE)),'my.model.')[[1]][2]
			dfMod[[cv]]$names1 <- my.model.RFRFE  
			names(dfMod[[cv]])[mod.ind[cv]] <- names1[1]
			mod.ind[cv] <- mod.ind[cv] +1 # update mod.ind indicator variable
		} # end CV types
    } # end RF-REF

    # --------------------------------------------
    # 8.11. SVM-RFE
    # --------------------------------------------
    if (fSVMRFE==T) {   # if SVM-RFE was selected, run the method
		outFile.SVMRFE <- file.path(PathDataSet,svmrfeFile)   # the same folder as the input is used for the output
      
		cat("-> SVM-RFE : Support Vector Machines Recursive Feature Elimination ...\n")
		# For each type of CV do all the statistics
		# -----------------------------------------------------
		for (cv in 1:length(CVtypes2)) {
		    cat("    -->",CVtypes2[cv],"\n")
		    ptmSVMRFE <- proc.time()
		    svmrfe.model <- SVMRFEreg(ds.train,ds.test,CVtypes[cv],i,fDet,outFile.SVMRFE,rfe_SVM_param_c,rfe_SVM_param_eps) # run SVM RFEet for each CV and regr method
		    print(proc.time() - ptmSVMRFE) # print running time
		    
		    my.stats.SVMRFE <- svmrfe.model$stat.values
		    my.model.SVMRFE <- svmrfe.model$model 
		    
		    #-------------------------------------------------------
		    # Add output from SVM RFE to the list of results
		    #-------------------------------------------------------
		    # List of results for each splitting, CV type & regression method
		    dfRes = mapply(c, my.stats.SVMRFE, dfRes, SIMPLIFY=F)
		    # List of models for each splitting, CV type & regression method
		    names1 <- strsplit(deparse(quote(my.model.SVMRFE)),'my.model.')[[1]][2]
		    dfMod[[cv]]$names1 <- my.model.SVMRFE  
		    names(dfMod[[cv]])[mod.ind[cv]] <- names1[1]
		    mod.ind[cv] <- mod.ind[cv] +1 # update mod.ind indicator variable
		  } # end CV types  
    } # end SVM RFE
    
    
    # Kill parallel server
    # ------------------------
    if (noCores!=1){
      # ------------------------------------------
      # parallel for windows:
      # ------------------------------------------
      if (Sys.info()[['sysname']]=="Windows"){ # clean the memory!
        stopCluster(cl)
      } 
    }
    # ------------------------
    

    # END OF REGRESSION METHODS/FUNCTIONS

    # -----------------------------------------------------------------------
    # (8.final) Produce comparison plots amongst models 
    # -----------------------------------------------------------------------

    for(cv in 1:length(CVtypes)){
      dfMod.n<- dfMod[[cv]]# keep only models with the same number of resamples
      dfMod.ind<- unlist(lapply(dfMod[[cv]],findResamps.funct))
      dfMod.ind.d<- which(duplicated(dfMod.ind)=='TRUE')[1]
      dfMod.in<- which(dfMod.ind!=dfMod.ind[dfMod.ind.d])
      dfMod.flag<- 0 # flag to indicate that dfMod consists of models with different resamples  
      
      if(is.na(dfMod.ind.d)!= TRUE){if(length(dfMod.in)!=0){dfMod.n<- dfMod.n[-dfMod.in]}}
      else{dfMod.flag<- 1}
      
      if(CVtypes[cv]!='LOOCV' && length(dfMod.n)>=2 && dfMod.flag!=1){ 
        cat("-> Comparisons plots for multiple regression methods ...\n")

        resamps <- resamples(dfMod.n)#,modelNames=names(dfMod[[cv]]))  
        # calculate their differences in terms of R2 and RMSE values
        difValues <- diff(resamps)
        #summary(difValues)
        #plot different models in terms of R2 adn RMSE values in the training set 
        pdf(file=paste(PathDataSet,"/ModelsComp.","iSplits.",i,".pdf",sep=""))
        print(bwplot(resamps, layout = c(2, 1),main=paste('Resampling results on the training set',' (data split ',i,')',sep='')))
        dev.off()    

        #plot differences of models in terms of R2 adn RMSE values in the training set 
        pdf(file=paste(PathDataSet,"/DifModels.R2.","iSplits.",i,".pdf",sep=""))
        print(dotplot(difValues,metric='Rsquared',main=paste('Models` differences on the training set',' (data split ',i,')',sep='')))
        dev.off()
        
        pdf(file=paste(PathDataSet,"/DifModels.RMSE.","iSplits.",i,".pdf",sep=""))
        print(dotplot(difValues,metric='RMSE',main=paste('Models` differences on the training set',' (data split ',i,')',sep='')))

        dev.off()
      }
    }
   
  } # END SPLITTING
  
  #------------------------------------------------------------------------------
  # 9. Results for all splittings (not ordered)
  #-------------------------------------------------------------------------------
  cat("-> Results for all splitings ...\n")
  df.res <- data.frame(dfRes)
  # print(df.res) # print all results as data frame
  
  # Writing the statistics into output files: one with detailed splits, other with only averages
  # File names includin paths for the statistics outputs (only averages and split detailed; +averages [to be implemented])
  ResBySplitsF <- file.path(PathDataSet,ResBySplits)   # the main output file with statistics for each split
  write.csv(df.res, file = ResBySplitsF)    # write statistics data frame into a CSV output file
  
  #-------------------------------------------------------------------------------------
  # Averaged values of the results by each Regression Method & CV type
  #-------------------------------------------------------------------------------------
  cat("-> Averaged statistics ...\n")
  ResAvgsF <- file.path(PathDataSet,ResAvgs)           # the main output file with averaged statistics for each regression method
  #library(data.table)
  dt.res  <- data.table(df.res) # convert data frame into data table (for sorting abd averaging)
  
  # MEANS for each Regression Method & CV type
  #--------------------------------------------------------------------------------------------------------------
  # means for all CV types, not only 10CV
  dt.mean <- dt.res[,list(adjR2.tr.Avg=mean(adjR2.tr),RMSE.tr.Avg=mean(RMSE.tr),R2.tr.Avg=mean(R2.tr),
                          RMSEsd.tr.Avg=mean(RMSEsd.tr),R2sd.tr.Avg=mean(R2sd.tr),adjR2.ts.Avg=mean(adjR2.ts),
                          RMSE.ts.Avg=mean(RMSE.ts),R2.ts.Avg=mean(R2.ts),corP.ts.Avg=mean(corP.ts),
                          adjR2.both.Avg=mean(adjR2.both),RMSE.both.Avg=mean(RMSE.both),
                          R2.both.Avg=mean(R2.both),NoModelFeats.Avg=round(mean(NoModelFeats),1)),by="RegrMeth,CVtype"]
  
  dt.mean     <- dt.mean[dt.mean$CVtype=="repeatedcv",]   # keep only the 10CV results to be used to find the best model
  
  if (fR2rule==T) { # R2 rule for the best model (default)
	dt.mean.ord <- dt.mean[order(-rank(R2.ts.Avg))]     # descendent order the averages by R2.ts.Avg
  } else { # adjR2 rule for the best model
	dt.mean.ord <- dt.mean[order(-rank(adjR2.ts.Avg))]  # descendent order the averages by adjR2.ts.Avg
  }
  
  # Write averages descendent ordered by R2.ts.Avg / adjR2.ts.Avg
  #-------------------------------------------------------------------------------
  write.csv(data.frame(dt.mean.ord), file = ResAvgsF)    # write statistics data frame into a CSV output file
  
  #------------------------------------------------------------------------------
  # 10. Best model selection - detailed statistics
  #-------------------------------------------------------------------------------
  cat("-> Best model analysis ...\n")
  
  # Algorithm to verify similar R2 / adjR2 values
  # -> from the best ones (+/- 0.05 of R2 / adjR2), chose the one with less variables, after that the one with min RMSE
  
  best.dt  <- dt.mean.ord[1] # the best model (R2.ts / adjR2.ts) should be the first value in the descendent ordered results

  # best.reg <- paste(best.dt$RegrMeth,collapse="") # best regression method
  
  # Best model rule: R2 or adjR2 for ordering
  # +/- 0.05 R2 / adjR2ts --> min(RMSE)
  # best.adjR2.ts is used for R2 or adjR2 rule
  if (fR2rule==T) { # R2 rule for the best model (default)
	best.adjR2.ts <- as.numeric(data.frame(best.dt)[,10]) # best R2.ts avgs
  } else { # adjR2 rule for the best model
	best.adjR2.ts <- as.numeric(data.frame(best.dt)[,8])  # best adjR2.ts avgs
  }
  
  # best model with R2 or adjR2.ts +/- 0.05 and min of RMSE for Avgs
  if (fR2rule==T) { # R2 rule for the best model (default)
	best.dt  <- dt.mean.ord[R2.ts.Avg %between% c(best.adjR2.ts-0.05,best.adjR2.ts+0.05)][which.min(RMSE.ts.Avg)]
  } else { # adjR2 rule for the best model
	best.dt  <- dt.mean.ord[adjR2.ts.Avg %between% c(best.adjR2.ts-0.05,best.adjR2.ts+0.05)][which.min(RMSE.ts.Avg)]
  }
  
  best.reg <- paste(best.dt$RegrMeth,collapse="") # best regrression method
  cat("    -> Method:",best.reg,"\n")

  # best model non-averaged ? no. of features
  # -----------------------------------------------
  # best.method <- dt.res[CVtype == "repeatedcv"][RegrMeth == best.reg] # best modes corresponding with the avg best values
  # best.method.mean <- mean(as.numeric(data.frame(best.method)[,11])) # best adjR2.ts)
  # dt.res[CVtype == "repeatedcv"][RegrMeth == best.reg][NoModelFeats == min(NoModelFeats)][RMSE.ts == min(RMSE.ts)]
  
  #----------------------------------------------------------
  # 11. Best model detailed statistics 
  #----------------------------------------------------------
  # Write the best model statistics
  ResBestF <- file.path(PathDataSet,ResBest)
  write.table("Averaged values for all spits: ",file=ResBestF,append=T,sep=",",col.names=F,row.names=F,quote=F)
  # write.csv(data.frame(best.dt), file = ResBestF)    # write statistics data frame into a CSV output file
  write.table(data.frame(best.dt), file=ResBestF,append=T,sep=",",col.names=T,quote=F) # write statistics data frame into a CSV output file
  
  # Use the last split for dataset (ds.train & ds.test) ! (or chose other one?)
  
  # Run the caret function with the method from the best method, for one training-test split only
  # and append the details in the best model output file
  
  # ----------------------------------
  # Parallel support - registed cores
  # ----------------------------------
  if (noCores!=1){
    # ------------------------------------------
    # parallel for Linux or Mac:
    # ------------------------------------------
    if (Sys.info()[['sysname']]=="Linux" | Sys.info()[['sysname']]=="Darwin"){
      registerDoMC(cores = noCores) # CPU cores
    }
    # ------------------------------------------
    # parallel for windows:
    # ------------------------------------------
    if (Sys.info()[['sysname']]=="Windows"){
     cl<-makeCluster(noCores,outfile="") 
      registerDoSNOW(cl)
    }   
  }

  if (best.reg=="lm") {
    my.stats.reg  <- LMreg(ds.train,ds.test,"repeatedcv",i,T,ResBestF)$stat.values # run GLM for each CV and regr method
  }
  if (best.reg=="glmStepAIC") {
    my.stats.reg  <- GLMreg(ds.train,ds.test,"repeatedcv",i,T,ResBestF)$stat.values # run GLM for each CV and regr method
  }
  if (best.reg=="pls") {
    my.stats.reg  <- PLSreg(ds.train,ds.test,"repeatedcv",i,T,ResBestF)$stat.values # run PLS for each CV and regr method
  }
  if (best.reg=="lasso.RMSE") {
    my.stats.reg  <- LASSOreg(ds.train,ds.test,"repeatedcv",i,T,ResBestF)$stat.values # run LASSO for each CV and regr method
  }
  if (best.reg=="svmRadial") {  
    my.stats.reg  <- SVRMreg(ds.train,ds.test,"repeatedcv",i,T,ResBestF,rfe_SVM_param_c)$stat.values # run SVRM Radial for each CV and regr method
  }
  if (best.reg=="nnet") {  
    my.stats.reg  <- NNreg(ds.train,ds.test,"repeatedcv",i,T,ResBestF)$stat.values # run NNet for each CV and regr method
  } 
  if (best.reg=="rf") {  
    my.stats.reg  <- RFreg(ds.train,ds.test,"repeatedcv",i,T,ResBestF)$stat.values # run NNet for each CV and regr method
  } 
  if (best.reg=="svmRFE") {  
    my.stats.reg  <- SVMRFEreg(ds.train,ds.test,"repeatedcv",i,T,ResBestF,rfe_SVM_param_c,rfe_SVM_param_eps)$stat.values # run NNet for each CV and regr method
  } 
  if (best.reg=="glmnet") {  
    my.stats.reg  <- ENETreg(ds.train,ds.test,"repeatedcv",i,T,ResBestF)$stat.values # run ENET for each CV and regr method
  }
  if (best.reg=="rfRFE") {  
    my.stats.reg  <- RFRFEreg(ds.train,ds.test,"repeatedcv",i,T,ResBestF)$stat.values # run RF RFE for each CV and regr method
  }

#   if (best.reg=="rbfDDA") {  
#     my.stats.reg  <- RBF_DDAreg(ds.train,ds.test,"repeatedcv",negThrStep,i,T,ResBestF)$stat.values # run rbfDDA for each CV and regr method
#  }

  #--------------------------------------------------------------------------------------
  # 12. Test best model with test dataset + Y randomization
  #--------------------------------------------------------------------------------------
  # ratios Yrand R2 - Best model R2 / Best model R2
  R2Diff.Yrand = NA
  if (noYrand >  0) {
    R2Diff.Yrand <- Yrandom(
      ds,trainFrac,best.reg,my.stats.reg$R2.ts,noYrand,ResBestF,rfe_SVM_param_c,rfe_SVM_param_eps
    ) # mean value of ratio (deatails are printed to output file)
  }
  
  # (ex param: negThrStep) 

  # Kill parallel server
  # ------------------------
  if (noCores!=1){
    # ------------------------------------------
    # parallel for windows:
    # ------------------------------------------
    if (Sys.info()[['sysname']]=="Windows"){ # clean the memory!
      stopCluster(cl)
    } 
  }
  # Assessment of Applicability Domain (plot leverage) was included as details in each regression function
  
  # Print total execution time
  cat("\nRRegrs total execution time\n")
  print(proc.time() - ptmTot) # print running time

  #----------------------------
  # Indicate main result files
  #----------------------------
  cat("\nMAIN RESULT FILES\n")
  cat("======================\n")
  cat("Statistics for all data set splittings/methods/CV types:", ResBySplits,"\n")
  cat("Averages by method/CV type:",ResAvgsF,"\n")
  cat("Best model statistics:",ResBestF,"\n")
  cat("Best model plots:",paste(ResBestF,".repeatedcv.split",i,".pdf",sep=""),"\n")
  if (noYrand > 0) cat("Best model Y-randomization plot:",paste(ResBestF,".Yrand.Hist.pdf",sep=""),"\n")
  cat("\n* if you choose Details, additional CSV and PDF files will be create for each method.\n")

  return(list(BestMethod=best.reg,BestStats=my.stats.reg, Models=dfMod))
  # return a list with 3 items: the name of the best method, the statistics for the best model, the list with all the fitted models (including the best one)
}
