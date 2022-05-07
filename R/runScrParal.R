
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
  maxDimensionScr=0.8*434*10*10^4,
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
  basicInfo=dataInfo(data=data,qualifyRefTax=TRUE,
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
    taxon_to_be_sample<-(results$goodRefTaxaCandi[!(results$goodRefTaxaCandi%in%refTaxa)])
    num_to_be_sample<-(nRef-length(refTaxa))
    if (num_to_be_sample>=length(taxon_to_be_sample)) {
      num_to_be_sample<-length(taxon_to_be_sample)
      message("The number of candidate reference taxon is smaller than the number of taxon required in phase 1.
              The number of taxon was set to be ",num_to_be_sample)
    }

    refTaxa_extra=sample(taxon_to_be_sample,num_to_be_sample)
    refTaxa=c(refTaxa,refTaxa_extra)
    if (length(refTaxa)==0) {
        stop("No candidate reference taxon is available. Please try to lower the reference taxon boundary.")
    }
  } else if(length(refTaxa)>nRef){
    if(length(seed)>0){
      set.seed(as.numeric(seed))
    }
    refTaxa=sample(refTaxa,nRef)
  }

  results$refTaxa=refTaxa


  #
  ## run original data screen
  #
  gc(reset=TRUE)

  screen1=originDataScreen(data=data,testCovInd=testCovInd,
                           nRef=nRef,refTaxa=refTaxa,
                           paraJobs=paraJobs,
                           allFunc=allFunc,Mprefix=Mprefix,
                           covsPrefix=covsPrefix,
                           binPredInd=binPredInd,
                           sequentialRun=sequentialRun,
                           adjust_method=adjust_method,
                           seed=seed)

  maxMemUsedInMb=sum(gc()[,6])

  results$countOfSelecForAPred=screen1$countOfSelecForAPred
  results$estOfSelectForAPred<-screen1$estOfSelectForAPred
  results$testCovCountMat=screen1$testCovCountMat
  results$testEstMat<-screen1$testEstMat




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

