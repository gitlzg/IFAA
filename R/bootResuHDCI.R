bootResuHDCI=function(
  data,
  refTaxa,
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

  penal=runBootLassoHDCI(x=x,y=y,nPredics=nPredics,nTaxa=nTaxa,
                         refTaxaPosition=ii,bootLassoAlpha=bootLassoAlpha,
                         bootB=bootB,standardize=standardize,
                         seed=seed)
  rm(x,y)

  results$finalBetaEst=penal$beta
  results$CIvecLow=penal$betaCIlow
  results$CIvecUp=penal$betaCIhi

  results$finalBetaEst.LPR=penal$beta.LPR
  results$CIvecLow.LPR=penal$betaCIlow.LPR
  results$CIvecUp.LPR=penal$betaCIhi.LPR

  rm(penal)
  return(results)
}


