##' @export


bootResuHDCI=function(
  data,
  refTaxa,
  independence=T,
  identity=T,
  bootB,
  bootLassoAlpha,
  binPredInd,
  covsPrefix,
  Mprefix
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
  dataForEst=dataRecovTrans(data=data,independence=independence,
                            identity=identity,ref=refTaxa,
                            Mprefix=Mprefix,covsPrefix=covsPrefix,
                            binPredInd=binPredInd)
  x=as(dataForEst$xTildalong,"sparseMatrix")
  y=as(dataForEst$UtildaLong,"sparseVector")
  rm(dataForEst)

  penal=runBootLassoHDCI(x=x,y=y,nPredics=nPredics,nTaxa=nTaxa,
                         refTaxaPosition=ii,bootLassoAlpha=bootLassoAlpha,bootB=bootB)
  rm(x,y)

  results$finalBetaEst=penal$beta
  results$CIvecLow=penal$betaCIlow
  results$CIvecUp=penal$betaCIhi
  rm(penal)
  return(results)
}


