
##' @export



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
  standardize=F,
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
  lamMax=max(abs(Matrix::colSums(as.matrix(x)*as.vector(y))))/nObsAll
  print("Start cross validatoin with picasso")
  cvStartTime= proc.time()[3]
  
  if(length(lambda)!=1){
    for(i in 1:nCV){
      print(paste("Start the",i,"th set of picasso cross validatoin"))
      cvStartTimei = proc.time()[3]
      
      cvResul=cvPicasso(x=x,y=y,lamMax=lamMax,lambda=lambda,nfolds=nfolds,zeroSDCut=zeroSDCut,
                        lambda.min.ratio=lambda.min.ratio,nLam=nLam,
                        method=method,family=family,standardize=standardize,seed=seed,seedi=seedi)
      cvExeTimei= (proc.time()[3] - cvStartTimei)/60
      print(paste("The",i,"th set of cross validation with picasso is done and took",cvExeTimei,"minutes"))
    }
    
    print(paste("cross validation with picasso is successful at the",i,"th cross validation"))
    
    cvAllTime= (proc.time()[3] - cvStartTime)/60
    print(paste("Cross validation with picasso is done and took", cvAllTime, "minutes"))
    
    results$lamList=cvResul$lamList
    results$SSE=cvResul$SSE
    lambda=cvResul$optiLam
  }
  
  rm(cvResul)
  
  finalRun=picasso(X=x,Y=y,lambda=lambda,method=method,family=family,standardize=standardize)
  endT=proc.time()[3]
  rm(x,y)
  
  # convert to a sparse vector format from sparse matrix format
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
  
  # parameter being set to zero
  if(length(xWithNearZeroSd)>0) {
    print("beta that should be set to 0: ")
    print(xWithNearZeroSd)
  }
  rm(xWithNearZeroSd)
  
  # return results
  results$betaInt=beta
  results$betaNoInt=beta[-seq(1,length(beta),by=(nPredics+1))]
  rm(beta,nPredics)
  results$lambda=lambda
  rm(lambda)
  
  results$successCV=1
  results$nCV=nCV
  rm(nCV)
  return(results)
}
