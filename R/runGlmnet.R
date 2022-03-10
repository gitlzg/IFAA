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
  standardize=FALSE,
  intercept=TRUE,
  zeroSDCut=10^(-20)
){

  results=list()


  nBeta=ncol(x)

  nObsAll=length(y)


  # remove near constant x columns
  sdX=apply(x,2,sd)
  xWithNearZeroSd=which(sdX<=zeroSDCut)
  if(length(xWithNearZeroSd)>0){
    x=x[,-xWithNearZeroSd,drop=FALSE]
  }
  rm(sdX)

  # calculate lambda max
  if(family=="gaussian"){lamMax=max(abs(Matrix::colSums(x*y)))/nObsAll}
  lamVec=seq(lamMax,0,length=(nLam+1))[1:nLam]

  cvStartTime= proc.time()[3]


    cvStartTimei = proc.time()[3]

    cvResul=cv.glmnet(x=x,y=as.vector(y),lambda=lamVec,nfolds=nfolds,
                      family=family,intercept=intercept,standardize=standardize)

    cvExeTimei= (proc.time()[3] - cvStartTimei)/60

    lamOpi=as.numeric(cvResul$lambda.min)
    cvm=as.vector(cvResul$cvm)*nObsAll
    nLamUsed=length(as.vector(cvResul$lambda))

    rm(cvResul)

  cvAllTime= (proc.time()[3] - cvStartTime)/60

  rm(lamVec)

  finalLassoRun=glmnet(x=x,y=as.vector(y),lambda=lamOpi,family=family,
                       intercept=intercept,standardize=standardize)
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
  rm(cvm)
  return(results)
}

