#======================================================================================================================
# General functions
#======================================================================================================================
mae <- function(pred,obs,na.rm=FALSE){mean(abs(pred - obs),na.rm=na.rm)}

r2.adj.t.funct<- function(obs,pred,num.pred){
  #obs==y, pred=predicted, num.pred=number of idependent variables (predictors)
  #t: traditional formula
  y.mean<- mean(obs)
  x.in<- sum((obs-pred)^2)/sum((obs-y.mean)^2)
  x.in<- 1-x.in #r squared
  
  x.in<- (1-x.in)*((length(obs)-1)/(length(obs)-num.pred-1))
  x.in<- 1 - x.in 
  return(x.in)
}
#----------------------------------------------------------------------------------------------------------------------
r2.adj.funct<- function(obs,pred,num.pred){
  #obs==y, pred=predicted, num.pred=number of idependent variables (predictors)
  
  x.in<- cor(obs,pred)^2
  x.in<- (1-x.in)*((length(obs)-1)/(length(obs)-num.pred-1))
  x.in<- 1 - x.in 
  return(x.in)
}
#----------------------------------------------------------------------------------------------------------------------
rmse.funct<- function(obs,pred){ 
  #obs==y, pred=predicted
  return(sqrt(mean((pred - obs)^2)))
}
#----------------------------------------------------------------------------------------------------------------------
r2.funct<- function(obs,pred){  
  #obs==y, pred=predicted
  x.in<- cor(obs,pred)^2
  return(x.in)
}
#----------------------------------------------------------------------------------------------------------------------
r2.t.funct<- function(obs,pred){  
  #obs==y, pred=predicted
  y.mean<- mean(obs)
  x.in<- sum((obs-pred)^2)/sum((obs-y.mean)^2)
  x.in<- 1-x.in #r squared
  return(x.in)
}
#----------------------------------------------------------------------------------------------------------------------