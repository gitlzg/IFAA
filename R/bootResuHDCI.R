#---------------------------------------------------------------------------------------
## bootstrap results function for Lasso OLS HDCI
#---------------------------------------------------------------------------------------

bootResuHDCI=function(
  data,
  refTaxa,
  originRefTaxNam,
  maxDimension=434*5*10^4,
  bootB,
  bootLassoAlpha,
  binPredInd,
  covsPrefix,
  Mprefix,
  paraJobs,
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
  
  maxSubSamplSiz=floor(maxDimension/xCol)
  nToSamplFrom=length(y)
  subSamplK=ceiling(nToSamplFrom/maxSubSamplSiz)
  if(subSamplK==1)maxSubSamplSiz=nToSamplFrom

  nRuns=(ceiling(subSamplK/3))
  
  availCores=max(1,(availableCores()-2))
  if(length(paraJobs)==0)paraJobs=availCores
  
  message(paraJobs, " parallel jobs are registered for bootstrapping in Phase 2.")
  
  for(k in 1:nRuns){
    rowToKeep=sample(nToSamplFrom,maxSubSamplSiz)
    xSub=as((x[rowToKeep,]),"sparseMatrix")
    ySub=as((y[rowToKeep]),"sparseVector")

    penal=runBootLassoHDCI(x=xSub,y=ySub,nPredics=nPredics,nTaxa=nTaxa,
                           refTaxaPosition=ii,bootLassoAlpha=bootLassoAlpha,
                           bootB=bootB,standardize=standardize,
                           paraJobs=paraJobs,
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
    if(k>0 & (k%%(ceiling(nRuns/10))==0)){
      message(floor(100*k/nRuns)," percent of the estimation analsis for the final reference taxon ",originRefTaxNam,
              " have been done.")
    }
  }
  rm(x,y)
  gc()
  
  results$finalBetaEst=finalBetaEst/nRuns
  results$CIvecLow=CIvecLow/nRuns
  results$CIvecUp=CIvecUp/nRuns
  
  results$finalBetaEst.LPR=finalBetaEst.LPR/nRuns
  results$CIvecLow.LPR=CIvecLow.LPR/nRuns
  results$CIvecUp.LPR=CIvecUp.LPR/nRuns
  
  return(results)
}