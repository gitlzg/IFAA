
#---------------------------------------------------------------------------------------
## parallel screen function to select independent reference taxa
#---------------------------------------------------------------------------------------

runScrParal=function(
  data,
  testCovInd,
  testCovInOrder,
  testCovInNewNam,
  nRef,
  paraJobs,
  refTaxa,
  fwerRate,
  maxDimensionScr=0.8*434*10*10^4,
  standardize,
  sequentialRun,
  allFunc=allFunc,
  refReadsThresh,
  SDThresh,
  SDquantilThresh,
  balanceCut,
  Mprefix,
  covsPrefix,
  binPredInd,
  adjust_method,
  seed){
  
  results=list()
  
  # load data info
  basicInfo=dataInfo(data=data,qualifyRefTax=T,
                     refReadsThresh=refReadsThresh,
                     SDThresh=SDThresh,SDquantilThresh=SDquantilThresh,
                     balanceCut=balanceCut,Mprefix=Mprefix,
                     covsPrefix=covsPrefix,
                     binPredInd=binPredInd)
  
  taxaNames=basicInfo$taxaNames
  nTaxa=basicInfo$nTaxa
  nPredics=basicInfo$nPredics
  nSub=basicInfo$nSub
  predNames=basicInfo$predNames
  
  results$goodRefTaxaCandi=basicInfo$goodRefTaxaCandi
  rm(basicInfo)
  gc()
  
  nNorm=nTaxa-1
  nAlphaNoInt=nPredics*nNorm
  nAlphaSelec=nPredics*nTaxa
  
  # make reference taxa list
  if(length(refTaxa)<nRef){
    if(length(seed)>0){
      set.seed(as.numeric(seed))
    }
    refTaxa_extra=sample((taxaNames[!(taxaNames%in%refTaxa)]),(nRef-length(refTaxa)))
    refTaxa=c(refTaxa,refTaxa_extra)
    results$refTaxa=refTaxa
  }
  
  if(length(refTaxa)>=nRef){
    if(length(seed)>0){
      set.seed(as.numeric(seed))
    }
    refTaxa=sample(refTaxa,nRef)
    results$refTaxa=refTaxa
  }
  
  #
  ## run original data screen
  #
  gc(reset=TRUE)
  
  screen1=originDataScreen(data=data,testCovInd=testCovInd,
                           nRef=nRef,refTaxa=refTaxa,fwerRate=fwerRate,
                           paraJobs=paraJobs,
                           allFunc=allFunc,Mprefix=Mprefix,
                           covsPrefix=covsPrefix,
                           binPredInd=binPredInd,
                           standardize=standardize,
                           sequentialRun=sequentialRun,
                           adjust_method=adjust_method,
                           seed=seed)
  
  maxMemUsedInMb=sum(gc()[,6])
  message("maximum memory used after phase 1a: ",maxMemUsedInMb," Mb")
  
  results$countOfSelecForAPred=screen1$countOfSelecForAPred
  results$estOfSelectForAPred<-screen1$estOfSelectForAPred
  yTildLongList=screen1$yTildLongList
  results$testCovCountMat=screen1$testCovCountMat
  results$testEstMat<-screen1$testEstMat
  # Alist=screen1$Alist
  # Wlist=screen1$Wlist
  # twoList=screen1$twoList
  # lengthTwoList=screen1$lengthTwoList
  
  
  
  rm(screen1)
  gc()
  
  nTestCov=length(testCovInd)
  results$nTestCov=nTestCov
  results$nTaxa=nTaxa
  results$nPredics=nPredics
  
  results$taxaNames=taxaNames
  rm(taxaNames)
  return(results)
}

