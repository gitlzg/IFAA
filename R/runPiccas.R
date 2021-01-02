

runPicasso=function(
  x,
  y,
  nPredics,
  family="gaussian",
  nfolds=10,
  nCV=1,
  cvRestartPosition=0.9,
  zeroSDCut=0,
  lambda=NULL,
  lambda.min.ratio=0.05,
  nLam=100,
  method="mcp",
  standardize,
  permutY,
  allFunc=allFunc,
  seed,
  seedi
){

  results=list()

  # check if there are x variables with zero variance
  nObsAll=length(y)
  nBeta=ncol(x)

  if(permutY){
    # permute y
    y=sample(y,size=nObsAll)
  }
  sdX=apply(x,2,sd)
  xWithNearZeroSd=which(sdX<=zeroSDCut)

  # remove near constant columns
  if(length(xWithNearZeroSd)>0){
    x=x[,-xWithNearZeroSd]
  }
  rm(sdX)
  nearZeroSd=length(xWithNearZeroSd)

  # calculate lambda max
  lamMax=max(abs(Matrix::colSums(x*y)))/nObsAll

  if(length(lambda)!=1){
    for(i in 1:nCV){
      cvResul=cvPicasso(x=x,y=y,lamMax=lamMax,lambda=lambda,nfolds=nfolds,zeroSDCut=zeroSDCut,
                        lambda.min.ratio=lambda.min.ratio,nLam=nLam,
                        method=method,family=family,standardize=standardize,seed=seed,seedi=seedi)
    }

    results$lamList=cvResul$lamList
    results$SSE=cvResul$SSE
    lambda=cvResul$optiLam
  }

  rm(cvResul)

  finalRun=picasso(X=x,Y=y,lambda=lambda,method=method,family=family,standardize=standardize)
  endT=proc.time()[3]
  rm(x,y)

  # convert to a sparse vector format from sparse matrix format

  OverallIntercp=finalRun$intercept[1]
  finalRunBeta=as(finalRun$beta[,1],"sparseVector")
  if(length(xWithNearZeroSd)>0){
    betaTrans=groupBetaToFullBeta(nTaxa=nBeta,nPredics=1,
                                  unSelectList=sort(xWithNearZeroSd),newBetaNoInt=finalRunBeta)
    beta=betaTrans$finalBeta
    rm(betaTrans)
  } else {
    beta=finalRunBeta
  }
  lambda=finalRun$lambda

  rm(finalRun,finalRunBeta)

  rm(xWithNearZeroSd)

  # return results
  results$betaInt=beta
  results$betaNoInt=beta[-seq(1,length(beta),by=(nPredics+1))]
  rm(beta,nPredics)
  results$lambda=lambda
  results$overalIntercp=OverallIntercp
  rm(lambda)

  results$successCV=1
  results$nCV=nCV
  rm(nCV)
  return(results)
}
