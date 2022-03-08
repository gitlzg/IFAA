allUserFunc=function(){
  return(c("dataRecovTrans", "AIcalcu","metaData",
           "dataSparsCheck","dataInfo","groupBetaToFullBeta",
           "runGlmnet","runlinear"))
}

calculate_se<-function(boot_CI_i) {
  se_est<-(boot_CI_i[2]-boot_CI_i[1])/(2*1.96)
  return(se_est)
}
