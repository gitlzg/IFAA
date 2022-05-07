#-------------------------------------------------------------------
## function for cross validation using linear model
#-------------------------------------------------------------------

runlinear=function(
  x,
  y,
  nPredics,
  fwerRate=0.25,
  adjust_method="fdr"
){

  results=list()

  nBeta=ncol(x)

  nObsAll=length(y)



  lm_res<-lm(as.vector(y)~as.matrix(x)-1)

  full_name_coef<-names(lm_res$coefficients)
  valid_coef<-summary(lm_res)$coefficients
  bootResu<-matrix(nrow = length(full_name_coef),ncol = 4)
  rownames(bootResu)<-full_name_coef
  bootResu[rownames(bootResu)%in%rownames(valid_coef),]<-valid_coef



  p_value_est<-bootResu[,4]
  p_value_est_noint<-p_value_est[-seq(1,length(p_value_est),by=(nPredics+1))]
  p_value_est_noint_adj<-p.adjust(p_value_est_noint,adjust_method)
  p_value_est_noint_adj[is.na(p_value_est_noint_adj)]<-1

  coef_est<-abs(bootResu[,1])
  coef_est_noint<-coef_est[-seq(1,length(coef_est),by=(nPredics+1))]
  coef_est_noint[is.na(coef_est_noint)]<-max(coef_est_noint,na.rm = TRUE)


  # return
  results$betaNoInt=p_value_est_noint_adj<fwerRate
  results$betaInt=p_value_est
  results$coef_est_noint=coef_est_noint


  return(results)
}

