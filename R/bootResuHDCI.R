
bootResuHDCI=function(
  data,
  refTaxa,
  maxDimension=5*434*10^4,
  bootB,
  bootLassoAlpha,
  binPredInd,
  covsPrefix,
  Mprefix,
  standardize,
  seed
){
  
  results=list()
  # load data info
  basicInfo=dataInfo(data=data,binPredInd=binPredInd,
                     covsPrefix=covsPrefix,Mprefix=Mprefix)
  taxaNames=basicInfo$taxaNames
  ii=which(basicInfo$taxaNames%in%refTaxa)
  
  nTaxa=basicInfo$nTaxa
  nPredics=basicInfo$nPredics
  rm(basicInfo)
  gc()
  
  nNorm=nTaxa-1
  nAlphaNoInt=nPredics*nNorm
  nAlphaSelec=nPredics*nTaxa
  
  countOfSelec=rep(0,nAlphaSelec)
  resultsByRefTaxon=list()
  
  # inital Lasso OLS estimate
  dataForEst=dataRecovTrans(data=data,ref=refTaxa,
                            Mprefix=Mprefix,covsPrefix=covsPrefix)
  
  x=as(dataForEst$xTildalong,"sparseMatrix")
  y=as(dataForEst$UtildaLong,"sparseVector")
  rm(dataForEst)
  
  xCol=ncol(x)
  subSamplSiz=maxDimension/xCol
  nToSamplFrom=length(y)
  if(nToSamplFrom<=subSamplSiz){
    subSamplK=1
    subSamplSiz=nToSamplFrom
  }
  if(nToSamplFrom>subSamplSiz & nToSamplFrom<=(2*subSamplSiz))subSamplK=2
  if(nToSamplFrom>(2*subSamplSiz))subSamplK=3
  cat("subSamplK:",subSamplK,"\n")
  for(k in 1:subSamplK){
    rowToKeep=sample(nToSamplFrom,subSamplSiz)
    xSub=as((x[rowToKeep,]),"sparseMatrix")
    ySub=as((y[rowToKeep]),"sparseVector")
    penal=runBootLassoHDCI(x=xSub,y=ySub,nPredics=nPredics,nTaxa=nTaxa,
                           refTaxaPosition=ii,bootLassoAlpha=bootLassoAlpha,
                           bootB=bootB,standardize=standardize,
                           seed=seed)
    rm(xSub,ySub)
    gc()
    finalBetaEst.k=penal$beta
    CIvecLow.k=penal$betaCIlow
    CIvecUp.k=penal$betaCIhi
    
    finalBetaEst.LPR.k=penal$beta.LPR
    CIvecLow.LPR.k=penal$betaCIlow.LPR
    CIvecUp.LPR.k=penal$betaCIhi.LPR
    rm(penal)
    gc()
    
    if(k==1){
      finalBetaEst=finalBetaEst.k
      CIvecLow=CIvecLow.k
      CIvecUp=CIvecUp.k
      
      finalBetaEst.LPR=finalBetaEst.LPR.k
      CIvecLow.LPR=CIvecLow.LPR.k
      CIvecUp.LPR=CIvecUp.LPR.k
    }
    if(k>1){
      finalBetaEst=finalBetaEst+finalBetaEst.k
      CIvecLow=CIvecLow+CIvecLow.k
      CIvecUp=CIvecUp+CIvecUp.k
      
      finalBetaEst.LPR=finalBetaEst.LPR+finalBetaEst.LPR.k
      CIvecLow.LPR=CIvecLow.LPR+CIvecLow.LPR.k
      CIvecUp.LPR=CIvecUp.LPR+CIvecUp.LPR.k
    }
  }
  rm(x,y)
  gc()
  
  results$finalBetaEst=finalBetaEst/subSamplK
  results$CIvecLow=CIvecLow/subSamplK
  results$CIvecUp=CIvecUp/subSamplK
  
  results$finalBetaEst.LPR=finalBetaEst.LPR/subSamplK
  results$CIvecLow.LPR=CIvecLow.LPR/subSamplK
  results$CIvecUp.LPR=CIvecUp.LPR/subSamplK
  
  return(results)
}

