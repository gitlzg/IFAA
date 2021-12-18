allUserFunc=function(){
  return(c("dataRecovTrans", "AIcalcu","metaData",
           "dataSparsCheck","dataInfo",
           "cvPicasso","runPicasso","groupBetaToFullBeta",
           "runGlmnet","runlinear"))
}

calculate_se<-function(boot_CI_i) {
  se_est<-(boot_CI_i[2]-boot_CI_i[1])/(2*1.96)
  # p_value_i<-pnorm(abs(est_i/se_est))
  return(se_est)
}
