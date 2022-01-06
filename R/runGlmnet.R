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
  zeroSDCut=0
  #zeroSDCut=10^(-40)
){

  results=list()

  # print("start glmnet function")
  # print(paste("nrow of x: ",nrow(x)))
  # print(paste("ncol of x: ",ncol(x)))

  nBeta=ncol(x)

  nObsAll=length(y)
  # print(paste("length of y: ",length(y)))

  # center y if continuous
  # if(family=="gaussian"){y=y-mean(y)}

  # remove near constant x columns
  sdX=apply(x,2,sd)
  xWithNearZeroSd=which(sdX<=zeroSDCut)
  if(length(xWithNearZeroSd)>0){
    x=x[,-xWithNearZeroSd,drop=FALSE]
  }
  rm(sdX)

  # calculate lambda max
  if(family=="gaussian"){lamMax=max(abs(Matrix::colSums(x*y)))/nObsAll}
  #lamVec=seq(lamMax,(lamMax*lambda.min.ratio),length=nLam)
  lamVec=seq(lamMax,0,length=(nLam+1))[1:nLam]

  # print("Start cross validatoin with glmnet")
  cvStartTime= proc.time()[3]

  #set.seed(123)

    cvStartTimei = proc.time()[3]

    cvResul=cv.glmnet(x=x,y=as.vector(y),lambda=lamVec,nfolds=nfolds,
                      family=family,intercept=intercept,standardize=standardize)

    cvExeTimei= (proc.time()[3] - cvStartTimei)/60
    # print(paste("The cross validation with glmnet is done and took",cvExeTimei,"minutes"))

    lamOpi=as.numeric(cvResul$lambda.min)
    cvm=as.vector(cvResul$cvm)*nObsAll
    nLamUsed=length(as.vector(cvResul$lambda))

    rm(cvResul)

  cvAllTime= (proc.time()[3] - cvStartTime)/60
  # print(paste("Cross validation with glmnet is done and took", cvAllTime, "minutes"))

  results$lamList=lamVec
  rm(lamVec)
  results$lambda=lamOpi
  rm(lamOpi)

  finalLassoRun=glmnet(x=x,y=as.vector(y),lambda=results$lambda,family=family,
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

    # print("beta should be set to 0:")
    # print(xWithNearZeroSd)
    rm(xWithNearZeroSd)
    } else {
     beta=finalLassoRunBeta
     }
  rm(finalLassoRun)

  # return
  results$betaNoInt=beta[-seq(1,length(beta),by=(nPredics+1))]
  results$betaInt=beta
  rm(beta)
  results$SSE=cvm
  rm(cvm)
  return(results)
}

# runGlmnet(x=x,y=y,nPredics)
