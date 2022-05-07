#-------------------------------------------------------------------
## function for cross validation using glmnet package
#-------------------------------------------------------------------

runGlmnet=function(
  x,
  y,
  nPredics,
  family="gaussian",
  nfolds=10,
  lambda.min.ratio=0.05,
  nLam=100,
  intercept=TRUE,
  zeroSDCut=10^(-20)
){

  results=list()

  nBeta=ncol(x)

  # remove near constant x columns
  sdX=apply(x,2,sd)
  xWithNearZeroSd=which(sdX<=zeroSDCut)
  if(length(xWithNearZeroSd)>0){
    x=x[,-xWithNearZeroSd,drop=FALSE]
  }
  rm(sdX)

    cvResul=cv.glmnet(x=x,y=as.vector(y),nlambda=nLam,nfolds=nfolds,
                      family=family,intercept=intercept,standardize=FALSE)

    lamOpi=as.numeric(cvResul$lambda.min)

    rm(cvResul)

  finalLassoRun=glmnet(x=x,y=as.vector(y),lambda=lamOpi,family=family,
                       intercept=intercept,standardize=FALSE)
  rm(x,y)

  # convert to regular matrix format
  finalLassoRunBeta=as.vector(finalLassoRun$beta)

  # convert back to the full beta if there near constant x columns
  if(length(xWithNearZeroSd)>0){
    betaTrans=groupBetaToFullBeta(nTaxa=nBeta,nPredics=1,
                                  unSelectList=sort(xWithNearZeroSd),newBetaNoInt=finalLassoRunBeta)
    beta=betaTrans$finalBeta
    rm(betaTrans)

    rm(xWithNearZeroSd)
    } else {
     beta=finalLassoRunBeta
     }
  rm(finalLassoRun)

  results$betaNoInt=beta[-seq(1,length(beta),by=(nPredics+1))]
  rm(beta)
  return(results)
}

