##' @export


getScrResu=function(
  method,
  data,
  testCovInd,
  nPermu,
  doPermut=T,
  nRef,
  paraJobs,
  refTaxa=NULL,
  fwerRate,
  independence,
  identity,
  goodIndeCutPerc=0.33,
  allFunc=allFunc,
  refReadsThresh,
  SDThresh,
  SDquantilThresh,
  balanceCut,
  Mprefix,
  covsPrefix,
  binPredInd
){
  results=list()

  # run permutation
  scrParal=runScrParal(data=data,testCovInd=testCovInd,nRef=nRef,
                       nPermu=nPermu,doPermut=doPermut,paraJobs=paraJobs,
                       refTaxa=refTaxa,independence=independence,identity=identity,
                       method=method,allFunc=allFunc,
                       refReadsThresh=refReadsThresh,SDThresh=SDThresh,
                       SDquantilThresh=SDquantilThresh,balanceCut=balanceCut,
                       Mprefix=Mprefix,covsPrefix=covsPrefix,binPredInd=binPredInd)

  selecMat=scrParal$countOfSelecForAPred
  testCovCountMat=scrParal$testCovCountMat
  taxaNames=scrParal$taxaNames
  goodRefTaxaCandi=scrParal$goodRefTaxaCandi
  nTaxa=scrParal$nTaxa
  nPredics=scrParal$nPredics
  nTestCov=scrParal$nTestCov
  results$selecMat=selecMat
  binomPar=scrParal$binomPar
  maxVec=scrParal$maxVec
  results$maxVec=maxVec
  permutTestCovByTaxa=scrParal$permutTestCovByTaxa
  rm(scrParal)
  gc()

  #print("maxVec:")
  #print(maxVec)

  results$taxaVecLeastCount=(selecMat[1,]==min(selecMat[1,]))+0
  results$taxaWithLeastCount=selecMat[,(selecMat[1,]==min(selecMat[1,]))]

  taxaLeastCountQualified=results$taxaWithLeastCount[names(results$taxaWithLeastCount)%in%goodRefTaxaCandi]

  if(length(taxaLeastCountQualified)>0){
    nObsNonZero=apply(data[,names(taxaLeastCountQualified),drop=F]>0,2,sum)
    results$taxaWithLeastCountMostNon0=taxaLeastCountQualified[tail(which(nObsNonZero==max(nObsNonZero)),n=1)]
  }
  if(length(taxaLeastCountQualified)==0){
    nObsNonZero=apply(data[,names(results$taxaWithLeastCount),drop=F]>0,2,sum)
    results$taxaWithLeastCountMostNon0=results$taxaWithLeastCount[tail(which(nObsNonZero==max(nObsNonZero)),n=1)]
    cat("Final reference taxon was not qualified, taxon with least selection times was used","\n")
  }
  rm(data)

  results$refTaxonQualified=length(taxaLeastCountQualified)>0

  goodRefWithCount=selecMat[,(colnames(selecMat)%in%goodRefTaxaCandi)]
  results$goodRefLeastCount=names(goodRefWithCount)[tail(which(goodRefWithCount==min(goodRefWithCount)),n=1)]

  if(length(refTaxa)>0){nRef=length(refTaxa)}

  if(doPermut){
    # control family wise error rate
    fwerCut=quantile(maxVec,probs=(1-fwerRate))
    results$selecTaxaFWER=(selecMat[1,]>=fwerCut)+0
    rm(maxVec)

    if(nTestCov==1)results$selecMatAll=matrix(results$selecTaxaFWER,nrow=1)

    if(nTestCov>1){
      fwerCutForAll=matrix(NA,nrow=nTestCov,ncol=nTaxa)
      for(i in 1:nTestCov){
        for(j in 1:nTaxa){
          fwerCutForAll[i,j]=quantile(permutTestCovByTaxa[[j]][i,],probs=(1-fwerRate))
        }
      }

      results$selecMatAll=(testCovCountMat>=fwerCutForAll)+0
      results$testCovCountMat=testCovCountMat
      results$fwerCutForAll=fwerCutForAll
      rm(permutTestCovByTaxa,testCovCountMat,fwerCutForAll)
    }

    goodIndpCut=quantile(selecMat[1,(selecMat[1,]<=fwerCut)],prob=goodIndeCutPerc)
    taxaLessGoodCut=selecMat[1,(selecMat[1,]<=goodIndpCut)]
    results$taxaNames=taxaNames

    rm(taxaNames)

    goodIndpRefTaxWithCount=taxaLessGoodCut[(names(taxaLessGoodCut)%in%goodRefTaxaCandi)]
    results$goodIndpRefTaxLeastCount=names(tail(goodIndpRefTaxWithCount[goodIndpRefTaxWithCount==min(goodIndpRefTaxWithCount)],n=1))

    rm(goodRefTaxaCandi,taxaLessGoodCut)
    results$fwerCut=fwerCut
  }

  if(!doPermut){results$goodIndpRefTaxLeastCount=refTaxa}
  return(results)
}
