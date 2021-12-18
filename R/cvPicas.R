
cvPicasso=function(
  x,
  y,
  family,
  lamMax,
  lambda,
  nfolds,
  zeroSDCut,
  lambda.min.ratio,
  nLam,
  method,
  standardize,
  seed,
  seedi
){

  results=list()

  nObsAll=length(y)
  nBeta=ncol(x)

  # obtain lambda range
  if(length(lambda)==0){
    #lamMin=lambda.min.ratio*lamMax
    lamMin=0
    lamList=seq(lamMax,lamMin,length=(nLam+1))[1:nLam]
  }
  if(length(lambda)>0){
    lamList=lambda
  }

  # partition the data randomly into nfolds partitions
  parSize=floor(nObsAll/nfolds)

  # if(length(seed)>0)set.seed(as.numeric(seed)+10^7+seedi)

  randomShuf=sample(nObsAll, nObsAll)
  sampleInd=list()
  for (i in 1:nfolds){
    if(i<nfolds){
      sampleInd[[i]]=randomShuf[(1+(i-1)*parSize):(i*parSize)]
    }else{
      sampleInd[[i]]=randomShuf[(1+(i-1)*parSize):nObsAll]
    }
  }
  rm(randomShuf)

  # cross validation
  cvPara=matrix(NA,nrow=nLam,ncol=nfolds)

  for(i in 1:nfolds){
    # check if there is zero-variance x in the partitions
    startT.i=proc.time()[3]
    xi=x[-sampleInd[[i]],]
    nObs.i=nrow(xi)
    sdX.i=apply(xi,2,sd)
    xWithNearZeroSd.i=which(sdX.i<=zeroSDCut)
    rm(sdX.i)

    # remove near constant columns in x
    if(length(xWithNearZeroSd.i)>0) {
      xi=xi[,-xWithNearZeroSd.i]
    }
    nearZeroSd.i=length(xWithNearZeroSd.i)

    yi=as(y[-sampleInd[[i]]],"sparseVector")

    cv.i=picasso(X=xi,Y=yi,lambda=lamList,method=method,family=family,standardize=standardize)
    rm(xi,yi)

    if(length(xWithNearZeroSd.i)>0) {
      missPositions=xWithNearZeroSd.i
      for(iLam in 2:nLam){
        missPositions=c(missPositions,xWithNearZeroSd.i+(iLam-1)*nBeta)
      }
      rm(xWithNearZeroSd.i)
      gc()
      betaa=as.vector(cv.i$beta)
      betaTrans=groupBetaToFullBeta(nTaxa=(nBeta*nLam),nPredics=1,
                                    unSelectList=sort(missPositions),newBetaNoInt=betaa)
      rm(betaa)
      gc()
      betaMatrix=Matrix(betaTrans$finalBeta,nrow=nBeta,ncol=nLam)
      rm(betaTrans)
    } else {betaMatrix=cv.i$beta}

    # training is done, start testing
    if(sum(is.na(cv.i$intercept))>0){
      cv.i$intercept[is.na(cv.i$intercept)]=mean(cv.i$intercept[!is.na(cv.i$intercept)])
    }
    if(sum(is.na(cv.i$intercept))>nLam*0.15){
      warnings("Many intercept estimates are missing in picasso estimates")
    }
    intcep=rep(1,length(sampleInd[[i]]))%*%t(as.matrix(cv.i$intercept))
    rm(cv.i)
    yHatTest.i=x[sampleInd[[i]],]%*%betaMatrix+intcep
    rm(betaMatrix,intcep)
    gc()

    resiVecs.i=(yHatTest.i-y[sampleInd[[i]]])
    rm(yHatTest.i)
    cvPara[,i]=apply(resiVecs.i,2,function(x)sum(x^2))
    rm(resiVecs.i)
  }
  rm(sampleInd,x,y)

  SSE=Matrix::rowSums(cvPara)
  rm(cvPara)

  optiLamInd=tail(which(SSE==min(SSE)),n=1)
  optiLam=lamList[optiLamInd]

  # return results
  results$SSE=SSE
  rm(SSE)

  results$lamList=lamList
  rm(lamList)

  results$optiLam=optiLam
  rm(optiLam)

  results$optiLamInd=optiLamInd
  rm(optiLamInd)

  results$nLam=nLam
  rm(nLam)
  return(results)
}

