
getScrResu=function(
  method,
  data,
  testCovInd,
  testCovInOrder,
  testCovInNewNam,
  nPermu,
  x1permut,
  nRef,
  paraJobs,
  refTaxa,
  fwerRate,
  goodIndeCutPerc=0.33,
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
  seed
){
  results=list()

  # run permutation
  scrParal=runScrParal(data=data,testCovInd=testCovInd,
                       testCovInOrder=testCovInOrder,
                       testCovInNewNam=testCovInNewNam,nRef=nRef,
                       nPermu=nPermu,
                       x1permut=x1permut,paraJobs=paraJobs,
                       refTaxa=refTaxa,method=method,
                       standardize=standardize,
                       sequentialRun=sequentialRun,
                       allFunc=allFunc,
                       refReadsThresh=refReadsThresh,SDThresh=SDThresh,
                       SDquantilThresh=SDquantilThresh,balanceCut=balanceCut,
                       Mprefix=Mprefix,covsPrefix=covsPrefix,
                       binPredInd=binPredInd,seed=seed)

  selecCountOverall=scrParal$countOfSelecForAPred
  selecCountMatIndv=scrParal$testCovCountMat
  taxaNames=scrParal$taxaNames
  goodRefTaxaCandi=scrParal$goodRefTaxaCandi

  nTaxa=scrParal$nTaxa
  nPredics=scrParal$nPredics
  nTestCov=scrParal$nTestCov
  binomPar=scrParal$binomPar
  maxVec=scrParal$maxVec
  results$maxVec=maxVec
  MaxMatTestCovByPermu=scrParal$MaxMatTestCovByPermu
  results$refTaxa=scrParal$refTaxa
  rm(scrParal)
  gc()

    # control family wise error rate
    originFwerCut=quantile(maxVec,probs=(1-fwerRate))
    fwerCut=originFwerCut

    results$twoMeanUsed=0

    if(fwerCut<=min(selecCountOverall)){
      # use two mean clustering
      twoMeanClusResu=twoMeanClus(matrix=selecCountOverall,nRef=nRef)
      fwerCut=twoMeanClusResu$twoMeanCut
      rm(twoMeanClusResu)
      if(fwerCut==0)fwerCut=0.1
      results$twoMeanUsed=1
    }
    results$fwerCut=fwerCut

    results$selecTaxaFWER=as(as.vector((selecCountOverall>=fwerCut)+0),"sparseVector")
    rm(maxVec)

    if(nTestCov==1){
      results$selecMatIndv=matrix(results$selecTaxaFWER,nrow=1)
      results$selecIndvInOverall=matrix(results$selecTaxaFWER,nrow=1)
      results$selecCountMatIndv=selecCountOverall
      results$fwerCutIndv=fwerCut
    }

    if(nTestCov>1){
      fwerCutIndv=vector()
      for(i in 1:nTestCov){
        originFwerCut.i=quantile(MaxMatTestCovByPermu[i,],probs=(1-fwerRate))
        fwerCutIndv[i]=originFwerCut.i
        if(originFwerCut.i<=min(selecCountMatIndv[i,])){
          # use two mean clustering
          twoMean=twoMeanClus(matrix=selecCountMatIndv[i,,drop=FALSE],nRef=nRef)
          fwerCutIndv[i]=twoMean$twoMeanCut
          if(twoMean$twoMeanCut==0)fwerCutIndv[i]=0.1
          rm(twoMean,originFwerCut.i)
        }
      }
      results$selecMatIndv=as((selecCountMatIndv>=fwerCutIndv)+0,"sparseMatrix")
      selecIndvInOverall_trans=t(as.matrix(results$selecMatIndv))*as.vector(results$selecTaxaFWER)
      results$selecIndvInOverall=as(t(as.matrix(selecIndvInOverall_trans)),"sparseMatrix")
      rm(selecIndvInOverall_trans)

      results$selecCountMatIndv=as(selecCountMatIndv,"sparseMatrix")
      results$fwerCutIndv=as(fwerCutIndv,"sparseVector")
      results$MaxMatTestCovByPermu=as(MaxMatTestCovByPermu,"sparseMatrix")
      rm(selecCountMatIndv,fwerCutIndv,MaxMatTestCovByPermu)
    }

    goodIndpCut=quantile(selecCountOverall[1,(selecCountOverall[1,]<=fwerCut)],prob=goodIndeCutPerc)
    taxaLessGoodCut=selecCountOverall[1,(selecCountOverall[1,]<=goodIndpCut)]
    taxaLessFWERCut=selecCountOverall[1,(selecCountOverall[1,]<fwerCut)]
    results$taxaNames=taxaNames
    results$selecCountOverall=selecCountOverall
    rm(taxaNames,selecCountOverall)

    goodIndpRefTaxWithCount=taxaLessGoodCut[(names(taxaLessGoodCut)%in%goodRefTaxaCandi)]
    if(length(goodIndpRefTaxWithCount)==0){
      results$goodIndpRefTaxLeastCount=NULL
    } else {
      results$goodIndpRefTaxLeastCount=names(tail(goodIndpRefTaxWithCount[goodIndpRefTaxWithCount==min(goodIndpRefTaxWithCount)],n=1))
    }

    results$goodIndpRefTaxWithCount=goodIndpRefTaxWithCount
    goodIndpRefTaxFWERcut=taxaLessFWERCut[(names(taxaLessFWERCut)%in%goodRefTaxaCandi)]
    results$goodIndpRefTaxFWERcut=goodIndpRefTaxFWERcut
    if(length(goodIndpRefTaxFWERcut)==0){
      results$goodIndpRefTaxFWERcutLeastCount=NULL
    } else {
      results$goodIndpRefTaxFWERcutLeastCount=names(tail(goodIndpRefTaxFWERcut[goodIndpRefTaxFWERcut==min(goodIndpRefTaxFWERcut)],n=1))
    }
    results$goodIndpRefTaxFWERcut=goodIndpRefTaxFWERcut
    results$goodRefTaxaCandi=goodRefTaxaCandi
    rm(goodRefTaxaCandi)

    minPosGoodCut=NULL
    if(sum(taxaLessGoodCut>0)>0){
      minPosGoodCut=min(taxaLessGoodCut[taxaLessGoodCut>0])
    }

    minPosFWERCut=NULL
    if(sum(taxaLessFWERCut>0)>0){
      minPosFWERCut=min(taxaLessFWERCut[taxaLessFWERCut>0])
    }
    results$taxaLessGoodCut=taxaLessGoodCut
    results$taxaLessFWERCut=taxaLessFWERCut

    results$minPosGoodCut=minPosGoodCut
    results$minPosFWERCut=minPosFWERCut

    results$refTaxonQualified=1
    if(length(results$goodIndpRefTaxLeastCount)==1){
      results$finalIndpRefTax=results$goodIndpRefTaxLeastCount
    } else if (length(results$goodIndpRefTaxFWERcutLeastCount)==1){
      results$finalIndpRefTax=results$goodIndpRefTaxFWERcutLeastCount
    } else if(length(minPosGoodCut)>0){
      results$finalIndpRefTax=names(tail(taxaLessGoodCut[taxaLessGoodCut==minPosGoodCut],n=1))
      results$refTaxonQualified=0
    } else if(length(minPosFWERCut)>0){
      results$finalIndpRefTax=names(tail(taxaLessFWERCut[taxaLessFWERCut==minPosFWERCut],n=1))
      results$refTaxonQualified=0
    }else{
      taxaLessGoodCutLeastCount=taxaLessGoodCut[taxaLessGoodCut==min(taxaLessGoodCut)]
      results$finalIndpRefTax=names(taxaLessGoodCutLeastCount)[1]
      results$refTaxonQualified=0
    }
    rm(taxaLessGoodCut,taxaLessFWERCut)
  
  return(results)
}
