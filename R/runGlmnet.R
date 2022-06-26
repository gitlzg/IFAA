#-------------------------------------------------------------------
## function for cross validation using glmnet package
#-------------------------------------------------------------------

runGlmnet = function(x,
                     y,
                     nPredics,
                     family = "gaussian",
                     intercept = TRUE) {
  results = list()
  
  ls.est = glmnet(
    x = x,
    y = as.vector(y),
    family = family,
    intercept = intercept,
    lambda=0,
    standardize = FALSE
    )
  
  finalLam=quantile(abs(ls.est$beta[ls.est$beta!=0]),probs=0.35)
  
  rm(ls.est)
  
  final.est = glmnet(
    x = x,
    y = as.vector(y),
    family = family,
    intercept = intercept,
    lambda=finalLam,
    standardize = FALSE
    )
  
  rm(x, y)
  
  beta <-as.vector(final.est$beta)
  
  rm(cvResul)

  results$betaNoInt = beta[-seq(1, length(beta), by = (nPredics + 1))]
  rm(beta)
  return(results)
}
