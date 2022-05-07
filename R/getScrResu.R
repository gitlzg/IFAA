
getScrResu=function(
  data,
  testCovInd,
  testCovInOrder,
  testCovInNewNam,
  nRef,
  paraJobs,
  refTaxa,
  goodIndeCutPerc=0.33,
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
  seed
){
  results=list()

  # run permutation
  scrParal=runScrParal(data=data,testCovInd=testCovInd,
                       testCovInOrder=testCovInOrder,
                       testCovInNewNam=testCovInNewNam,nRef=nRef,
                       paraJobs=paraJobs,
                       refTaxa=refTaxa,
                       sequentialRun=sequentialRun,
                       allFunc=allFunc,
                       refReadsThresh=refReadsThresh,SDThresh=SDThresh,
                       SDquantilThresh=SDquantilThresh,balanceCut=balanceCut,
                       Mprefix=Mprefix,covsPrefix=covsPrefix,
                       binPredInd=binPredInd,
                       adjust_method=adjust_method,seed=seed)

  selecCountOverall=scrParal$countOfSelecForAPred
  selecEstOverall<-scrParal$estOfSelectForAPred

  selecCountMatIndv=scrParal$testCovCountMat
  selecEstMatIndv=scrParal$testEstMat

  taxaNames=scrParal$taxaNames
  goodRefTaxaCandi=scrParal$goodRefTaxaCandi

  nTaxa=scrParal$nTaxa
  nPredics=scrParal$nPredics
  nTestCov=scrParal$nTestCov
  results$refTaxa=scrParal$refTaxa
  rm(scrParal)
  gc()

  if(nTestCov==1){
    results$selecCountMatIndv=selecCountOverall
    results$selecEstMatIndv=selecEstOverall
  }

  if(nTestCov>1){
    results$selecCountMatIndv=as(selecCountMatIndv,"sparseMatrix")
    results$selecEstMatIndv=as(selecEstMatIndv,"sparseMatrix")

    rm(selecCountMatIndv)
  }
  goodIndpRefTaxWithCount=selecCountOverall[1,(colnames(selecCountOverall)%in%goodRefTaxaCandi)]
  goodIndpRefTaxWithEst=selecEstOverall[1,(colnames(selecEstOverall)%in%goodRefTaxaCandi)]
  if(length(goodIndpRefTaxWithCount)==0){
    results$goodIndpRefTaxLeastCount=NULL
  } else {
    results$goodIndpRefTaxLeastCount<-
      names(goodIndpRefTaxWithCount)[order(goodIndpRefTaxWithCount,abs(goodIndpRefTaxWithEst))][1:2]
    goodIndpRefTaxWithEst<-abs(goodIndpRefTaxWithEst[order(goodIndpRefTaxWithCount,abs(goodIndpRefTaxWithEst))])
    goodIndpRefTaxWithCount<-goodIndpRefTaxWithCount[order(goodIndpRefTaxWithCount,abs(goodIndpRefTaxWithEst))]
  }

  results$selecCountOverall=selecCountOverall

  results$goodIndpRefTaxWithCount=goodIndpRefTaxWithCount
  results$goodIndpRefTaxWithEst<-goodIndpRefTaxWithEst

  results$goodRefTaxaCandi=goodRefTaxaCandi
  rm(goodRefTaxaCandi)

  results$refTaxonQualified=2
  results$finalIndpRefTax=results$goodIndpRefTaxLeastCount

  return(results)
}
