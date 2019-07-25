
##' @export


dataRecovTrans=function(
  data,
  ref,
  Mprefix,
  covsPrefix,
  binPredInd,
  independence,
  identity
){
  results=list()

  # load A and log-ratio transformed RA
  data.and.init=dataAndInit(data=data,ref=ref,Mprefix=Mprefix,
                            covsPrefix=covsPrefix,
                            binPredInd=binPredInd,
                            independence=independence,identity=identity)
  rm(data)

  taxaNames=data.and.init$taxaNames
  xtremeSparsTaxa=data.and.init$xtremeSparsTaxa
  A=data.and.init$A
  logRatiow=data.and.init$logRatiow
  nSub=data.and.init$nSub
  nTaxa=data.and.init$nTaxa
  xData=data.and.init$xData
  nPredics=data.and.init$nPredics
  twoList=data.and.init$twoList
  lengthTwoList=data.and.init$lengthTwoList
  cov=data.and.init$cov
  taxaLowSd=data.and.init$taxaLowSd

  if(length(binPredInd)>0){
    results$taxaAndBinIndex=data.and.init$taxaAndBinIndex
  }
  rm(data.and.init)

  nNorm=nTaxa-1
  xDimension=nPredics+1 # predictors+intercept

  # create omegaRoot
  omegaRoot=list()
  for (j in 1:lengthTwoList){
    i=twoList[j]
    omegaRoot[[i]]=sqrtm(solve(A[[i]]%*%cov%*%t(A[[i]])))
  }
  rm(cov)

  # create X_i in the regression equation using Kronecker product
  xDataWithInter=data.matrix(cbind(rep(1,nSub),xData))
  colnames(xDataWithInter)[1]="Inter"
  rm(xData)

  xInRegres=list()
  for (j in 1:lengthTwoList){
    i=twoList[j]
    xInRegres[[i]]=t(kronecker(diag(nNorm),xDataWithInter[i,]))
  }

  rm(xDataWithInter)

  # transfer U and xInRegres using omega root
  xDataTilda=list()
  Utilda=list()

  for (j in 1:lengthTwoList){
    i=twoList[j]
    Utilda[[i]]=omegaRoot[[i]]%*%logRatiow[[i]]
    xDataTilda[[i]]=omegaRoot[[i]]%*%A[[i]]%*%xInRegres[[i]]
  }
  rm(omegaRoot,logRatiow,xInRegres)

  # stack Utilda
  for (j in 1:lengthTwoList){
    i=twoList[j]
    if (j==1) {UtildaLong=Utilda[[i]]
    } else {
      UtildaLong=rbind(UtildaLong,Utilda[[i]])
    }
  }
  rm(Utilda)

  # stack xTilda
  for (j in 1:lengthTwoList){
    i=twoList[j]
    if (j==1) {xTildalong=xDataTilda[[i]]
    } else {
      xTildalong=rbind(xTildalong,xDataTilda[[i]])
    }
  }
  rm(xDataTilda,lengthTwoList)


  # set the extreme spars taxa having zero x so glmnet set the coefficients to be 0
  if(length(xtremeSparsTaxa)>0){
    xtremeTaxaNums=which(taxaNames%in%xtremeSparsTaxa)
    startPositions=(xtremeTaxaNums-1)*(nPredics+1)+2
    stopPositions=xtremeTaxaNums*(nPredics+1)
    for (i in 1:length(xtremeTaxaNums)){
      xTildalong[,(startPositions[i]:stopPositions[i])]=0
    }
    rm(xtremeTaxaNums,startPositions,stopPositions)
  }

  # set the low sd taxa having zero x so glmnet set the coefficients to be 0
  if(length(taxaLowSd)>0){
    lowSDTaxaNums=which(taxaNames%in%taxaLowSd)
    startPositions=(lowSDTaxaNums-1)*(nPredics+1)+2
    stopPositions=lowSDTaxaNums*(nPredics+1)
    for (i in 1:length(lowSDTaxaNums)){
      xTildalong[,(startPositions[i]:stopPositions[i])]=0
    }
    rm(lowSDTaxaNums,startPositions,stopPositions)
  }

  # if there are binary predicttors, set the taxaAndBin columns to be 0 for glmnet
  if(length(binPredInd)>0){
    xTildalong[,results$taxaAndBinIndex]=0
  }

  # return objects
  results$UtildaLong=UtildaLong
  rm(UtildaLong)
  results$xTildalong=xTildalong
  rm(xTildalong)

  results$taxaNames=taxaNames
  rm(taxaNames)
  results$xtremeSparsTaxa=xtremeSparsTaxa
  rm(xtremeSparsTaxa)
  results$taxaLowSd=taxaLowSd
  rm(taxaLowSd)
  results$binPredInd=binPredInd
  rm(binPredInd)
  return(results)
}
