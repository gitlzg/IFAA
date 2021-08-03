
Regulariz=function(
  data,
  testCovInd,
  testCovInOrder,
  testCovInNewNam,
  microbName,
  nRef,
  nRefMaxForEsti,
  nPermu,
  x1permut,
  refTaxa,
  paraJobs,
  binaryInd,
  covsPrefix,
  Mprefix,
  reguMethod,
  fwerRate,
  bootB,
  bootLassoAlpha,
  standardize,
  sequentialRun,
  allFunc=allFunc,
  refReadsThresh,
  SDThresh,
  SDquantilThresh,
  balanceCut,
  seed
){
  results=list()
  regul.start.time = proc.time()[3]

  nTestCov=length(testCovInd)

  dataSparsCheck(data=data,Mprefix=Mprefix)

  # load abundance data info
  data.info=dataInfo(data=data,Mprefix=Mprefix,
                     covsPrefix=covsPrefix,
                     binPredInd=binaryInd)
  nSub=data.info$nSub
  taxaNames=data.info$taxaNames
  nPredics=data.info$nPredics
  nTaxa=data.info$nTaxa
  rm(data.info)

  results$x1permut=x1permut

  regul.start.time = proc.time()[3]
  message("Start Phase 1 association identification")
  selectRegroup=getScrResu(data=data,testCovInd=testCovInd,
                           testCovInOrder=testCovInOrder,
                           testCovInNewNam=testCovInNewNam,nRef=nRef,
                           nPermu=nPermu,paraJobs=paraJobs,
                           x1permut=x1permut,
                           refTaxa=refTaxa,fwerRate=fwerRate,
                           method=reguMethod,
                           standardize=standardize,
                           sequentialRun=sequentialRun,
                           allFunc=allFunc,
                           refReadsThresh=refReadsThresh,
                           SDThresh=SDThresh,
                           SDquantilThresh=SDquantilThresh,
                           balanceCut=balanceCut,
                           Mprefix=Mprefix,
                           covsPrefix=covsPrefix,
                           binPredInd=binaryInd,seed=seed)

  results$selecTaxaFWER=selectRegroup$selecTaxaFWER
  results$selecCountOverall=selectRegroup$selecCountOverall
  results$selecMatIndv=selectRegroup$selecMatIndv
  results$selecIndvInOverall=selectRegroup$selecIndvInOverall
  results$selecCountMatIndv=selectRegroup$selecCountMatIndv
  results$maxVec=selectRegroup$maxVec
  results$MaxMatTestCovByPermu=selectRegroup$MaxMatTestCovByPermu

  finalIndpRefTax=microbName[taxaNames%in%(selectRegroup$finalIndpRefTax)]
  results$finalRefTaxonQualified=selectRegroup$refTaxonQualified

  results$goodIndpRefTaxLeastCount=microbName[taxaNames%in%(selectRegroup$goodIndpRefTaxLeastCount)]
  results$goodIndpRefTaxWithCount=selectRegroup$goodIndpRefTaxWithCount
  names(results$goodIndpRefTaxWithCount)=microbName[taxaNames%in%names(selectRegroup$goodIndpRefTaxWithCount)]

  results$goodIndpRefTaxFWERcutLeastCount=microbName[taxaNames%in%(selectRegroup$goodIndpRefTaxFWERcutLeastCount)]
  results$goodIndpRefTaxFWERcut=selectRegroup$goodIndpRefTaxFWERcut
  names(results$goodIndpRefTaxFWERcut)=microbName[taxaNames%in%names(selectRegroup$goodIndpRefTaxFWERcut)]

  results$goodRefTaxaCandi=microbName[taxaNames%in%(selectRegroup$goodRefTaxaCandi)]

  results$taxaLessGoodCut=selectRegroup$taxaLessGoodCut
  names(results$taxaLessGoodCut)=microbName[taxaNames%in%names(selectRegroup$taxaLessGoodCut)]


  results$taxaLessFWERCut=selectRegroup$taxaLessFWERCut
  names(results$taxaLessFWERCut)=microbName[taxaNames%in%names(selectRegroup$taxaLessFWERCut)]

  results$randomRefTaxa=microbName[taxaNames%in%(selectRegroup$refTaxa)]

  results$fwerCut=selectRegroup$fwerCut
  results$twoMeanUsed=selectRegroup$twoMeanUsed
  results$fwerCutIndv=selectRegroup$fwerCutIndv
  rm(selectRegroup)

  MCPExecuTime = (proc.time()[3] - regul.start.time)/60
  results$MCPExecuTime=MCPExecuTime
  message("Phase 1 Associaiton identification is done and used ", round(MCPExecuTime,2)," minutes")

  results$finalizedBootRefTaxon=finalIndpRefTax

  startT=proc.time()[3]
  message("Start Phase 2 parameter estimation")

  unestimableTaxa=c()
  qualifyData=data
  if(length(binaryInd)>0){
    qualifyData=data[rowSums(data[,results$taxaNames]>0)>=2,,drop=FALSE]
    firstBinPredNam=paste0(covsPrefix,binaryInd)
    allBinPred=paste0(covsPrefix,binaryInd:nPredics)
    nBinPred=length(allBinPred)

    # find the pairs of binary preds and taxa for which the assocaiton is not identifiable
    AllTaxaNamesNoRefTax=taxaNames[!taxaNames%in%(results$finalizedBootRefTaxon)]
    unbalanceTaxa=c()
    for(i in AllTaxaNamesNoRefTax){
      for(j in allBinPred){
        twoColumns.ij=qualifyData[,c(i,j)]
        nNonZero=sum(twoColumns.ij[,1]>0)
        sumOfBin=sum(twoColumns.ij[(twoColumns.ij[,1]>0),2])
        if(sumOfBin%in%c(0,1,(nNonZero-1),nNonZero)){
          unbalanceTaxa=c(unbalanceTaxa,i)
        }
      }
    }
    unestimableTaxa=unique(unbalanceTaxa)
    rm(allBinPred,unbalanceTaxa)
  }

  # check zero taxa and subjects with zero taxa reads
  TaxaNoReads=which(Matrix::colSums(qualifyData[,taxaNames])==0)
  rm(qualifyData)
  unestimableTaxa=unique(c(unestimableTaxa,taxaNames[TaxaNoReads]))
  results$unEstTaxaPos=which(taxaNames%in%unestimableTaxa)
  rm(TaxaNoReads,unestimableTaxa)

  goodIndpRefTax.ascend=sort(results$goodIndpRefTaxWithCount)

  goodIndpRefTaxNam=names(goodIndpRefTax.ascend)
  allRefTaxNam=unique(c(results$finalizedBootRefTaxon,goodIndpRefTaxNam))
  nGoodIndpRef=length(allRefTaxNam)
  results$allRefTaxNam=allRefTaxNam

  results$nRefUsedForEsti=min(nGoodIndpRef,nRefMaxForEsti)

  message("Final Reference Taxa are: ",allRefTaxNam[seq(results$nRefUsedForEsti)])

  results$estiList=list()
  for(iii in 1:(results$nRefUsedForEsti)){
    message("Start estimation for the ", iii,"th final reference taxon: ",allRefTaxNam[iii])
    time11=proc.time()[3]
    originTaxNam=allRefTaxNam[iii]
    newRefTaxNam=taxaNames[microbName%in%originTaxNam]
    results$estiList[[originTaxNam]]=bootResuHDCI(data=data,
                                                       refTaxa=newRefTaxNam,
                                                       originRefTaxNam=originTaxNam,
                                                       bootB=bootB,bootLassoAlpha=bootLassoAlpha,
                                                       binPredInd=binaryInd,covsPrefix=covsPrefix,
                                                       Mprefix=Mprefix,paraJobs=paraJobs,
                                                       standardize=standardize,
                                                       seed=seed)
    time12=proc.time()[3]
    message("Estimation done for the ", iii,"th final reference taxon: ",allRefTaxNam[iii],
        " and it took ",round((time12-time11)/60,3)," minutes")
  }
  estiResults=results$estiList[[results$finalizedBootRefTaxon]]

  endT=proc.time()[3]

  message("Phase 2 parameter estimation done and took ",round((endT-startT)/60,3)," minutes.")


  results$betaMat=as(matrix(estiResults$finalBetaEst,nrow=nPredics),"sparseMatrix")
  results$CILowMat=as(matrix(estiResults$CIvecLow,nrow=nPredics),"sparseMatrix")
  results$CIUpMat=as(matrix(estiResults$CIvecUp,nrow=nPredics),"sparseMatrix")

  results$betaMat.LPR=as(matrix(estiResults$finalBetaEst.LPR,nrow=nPredics),"sparseMatrix")
  results$CILowMat.LPR=as(matrix(estiResults$CIvecLow.LPR,nrow=nPredics),"sparseMatrix")
  results$CIUpMat.LPR=as(matrix(estiResults$CIvecUp.LPR,nrow=nPredics),"sparseMatrix")
  rm(estiResults)

  estByCovList=list()

  for(i in 1:nTestCov){
    sigTaxaPosition=which(results$selecIndvInOverall[i,]!=0)
    nrow=length(sigTaxaPosition)
    if(nrow==0)next

    ncol=3
    estByCovMat=matrix(NA,nrow=nrow,ncol=ncol)

    for(j in 1:nrow){
      for(k in 1:ncol){
        if(k==1)estByCovMat[j,k]=results$betaMat.LPR[i,sigTaxaPosition[j]]
        if(k==2)estByCovMat[j,k]=results$CILowMat.LPR[i,sigTaxaPosition[j]]
        if(k==3)estByCovMat[j,k]=results$CIUpMat.LPR[i,sigTaxaPosition[j]]
      }
    }

    rownames(estByCovMat)=microbName[results$selecIndvInOverall[i,]!=0]
    colnames(estByCovMat)=c("Beta.LPR","LowB95%CI.LPR","UpB95%CI.LPR")

    rowsToKeep=which(estByCovMat[,1]!=0)
    estByCovMat=estByCovMat[rowsToKeep, , drop = FALSE]

    estByCovList[[testCovInOrder[i]]]=estByCovMat
    rm(estByCovMat)
  }

  if(length(estByCovList)==0){
    results$estByCovList="No significant assoication is identified. You may consider a higher FDR level."
  }else{
    results$estByCovList=estByCovList
  }
  rm(estByCovList)

  SigCovByTaxaList=list()
  for(i in 1:nTaxa){
    sigCov=which(results$selecIndvInOverall[,i]!=0)
    if(length(sigCov)==0)next
    SigCovByTaxaList[[microbName[i]]]=testCovInOrder[sigCov]
  }
  if(length(SigCovByTaxaList)==0){
    results$SigCovByTaxaList="No significant assoication is identified. You may consider a higher FDR level."
  }else{
    results$SigCovByTaxaList=SigCovByTaxaList
  }
  rm(SigCovByTaxaList,microbName)

  results$reguMethod=reguMethod
  results$nSub=nSub
  results$nTaxa=nTaxa
  results$nPredics=nPredics
  return(results)

  rm(data)

  # return results
  results$fwerRate=fwerRate
  results$nPermu=nPermu
  results$nRef=nRef
  return(results)
}

#
