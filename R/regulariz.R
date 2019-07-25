
##' @export



Regulariz=function(
  data,
  testCovInd,
  testCovInOrder,
  microbName,
  nRef,
  nPermu,
  refTaxa,
  paraJobs,
  binaryInd,
  covsPrefix,
  Mprefix,
  independence,
  identity,
  reguMethod,
  fwerRate,
  boot,
  bootB,
  bootLassoAlpha,
  allFunc=allFunc,
  refReadsThresh,
  SDThresh,
  SDquantilThresh,
  balanceCut
){
  results=list()
  regul.start.time = proc.time()[3]

  nTestCov=length(testCovInd)

  dataSparsCheck(data=data,Mprefix=Mprefix)

  # load abundance data info
  dataForEst=dataInfo(data=data,Mprefix=Mprefix,
                      covsPrefix=covsPrefix,
                      binPredInd=binaryInd)
  nSub=dataForEst$nSub
  taxaNames=dataForEst$taxaNames
  EName=dataForEst$EName
  nPredics=dataForEst$nPredics
  nTaxa=dataForEst$nTaxa
  rm(dataForEst)

  # run MCP penalization
  regul.start.time = proc.time()[3]
  cat("Start association identification","\n")
  #if(length(refTaxa)==0){
  selectRegroup=getScrResu(data=data,testCovInd=testCovInd,nRef=nRef,
                           nPermu=nPermu,paraJobs=paraJobs,
                           refTaxa=NULL,fwerRate=fwerRate,
                           independence=independence,identity=identity,
                           method=reguMethod,allFunc=allFunc,
                           refReadsThresh=refReadsThresh,
                           SDThresh=SDThresh,
                           SDquantilThresh=SDquantilThresh,
                           balanceCut=balanceCut,
                           Mprefix=Mprefix,
                           covsPrefix=covsPrefix,
                           binPredInd=binaryInd)

  results$selecTaxaFWER=selectRegroup$selecTaxaFWER
  results$countMat=selectRegroup$selecMat
  results$selecMatAll=selectRegroup$selecMatAll
  results$testCovCountMat=selectRegroup$testCovCountMat
  if(nTestCov==1)results$testCovCountMat=results$countMat
  results$maxVec=selectRegroup$maxVec
  results$goodIndpRefTaxLeastCount=selectRegroup$goodIndpRefTaxLeastCount
  results$goodRefLeastCount=selectRegroup$goodRefLeastCount

  results$taxaVecLeastCount=selectRegroup$taxaVecLeastCount
  results$taxaWithLeastCount=selectRegroup$taxaWithLeastCount
  results$taxaWithLeastCountMostNon0=selectRegroup$taxaWithLeastCountMostNon0
  results$finalRefTaxonQualified=selectRegroup$refTaxonQualified

  results$fwerCut=selectRegroup$fwerCut
  results$fwerCutForAll=selectRegroup$fwerCutForAll
  results$taxaNames=selectRegroup$taxaNames
  rm(selectRegroup)

  MCPExecuTime = (proc.time()[3] - regul.start.time)/60
  results$MCPExecuTime=MCPExecuTime
  cat("Associaiton identification is done and used", MCPExecuTime,"minutes","\n")

  cat("Deciding the final reference taxon","\n")

  if(length(refTaxa)==0){
    results$truGoodIndpRefTax=(length(results$goodIndpRefTaxLeastCount)>0)+0

    results$goodIndpRefTaxLeastCount=names(results$taxaWithLeastCountMostNon0)

    if(length(results$goodIndpRefTaxLeastCount)==0){
      if(length(results$goodRefLeastCount)>0){
        warning("An independent reference taxon was not found. The good ref taxon with least count of association is used instead.")
        results$goodIndpRefTaxLeastCount=results$goodRefLeastCount
      }
      if(length(results$goodRefLeastCount)==0){
        warning("An independent reference or good reference taxon was not found. The taxon with least positive count of association is used instead.")
        minPosCount=min(results$selecMat[results$selecMat>0])
        results$goodIndpRefTaxLeastCount=results$taxaNames[tail(which(results$selecMat==minPosCount),n=1)]
      }
    }

    results$finalizedBootRefTaxon=results$goodIndpRefTaxLeastCount

    results$finalRefTaxon=microbName[results$taxaNames==results$finalizedBootRefTaxon]

    cat("Final Reference Taxon is",results$finalRefTaxon,"\n")

    if(results$finalRefTaxonQualified==0){
      cat("Reference taxon was not qualified, taxon with least selection times was used","\n")
    }

    cat("Start association estimation","\n")

    estiResults=bootResuHDCI(data=data,refTaxa=results$goodIndpRefTaxLeastCount,independence=independence,
                             identity=identity,
                             bootB=bootB,bootLassoAlpha=bootLassoAlpha,
                             binPredInd=binaryInd,
                             covsPrefix=covsPrefix,Mprefix=Mprefix)
    cat("Association estimation done","\n")
  }
  if(length(refTaxa)>0){
    CItime1 = proc.time()[3]
    refTaxaWithNewNam=paste0(Mprefix,which(microbName==refTaxa))
    results$finalRefTaxon=refTaxaWithNewNam

    cat("Final Reference Taxon is",results$finalRefTaxon,"\n")

    estiResults=list()
    cat("Start association estimation","\n")
    estiResults=bootResuHDCI(data=data,refTaxa=refTaxaWithNewNam,
                             independence=independence,
                             identity=identity,
                             bootB=bootB,bootLassoAlpha=bootLassoAlpha,
                             binPredInd=binaryInd,
                             covsPrefix=covsPrefix,Mprefix=Mprefix)

    cat("Association estimation done","\n")
    CItime2 = proc.time()[3]
    CIduration=(CItime2-CItime1)/60
    print(paste("Calculation of CI took", CIduration,"minutes"))
  }

  results$estiResults=estiResults
  results$betaMat=as(matrix(estiResults$finalBetaEst,nrow=nPredics),"sparseMatrix")
  results$CILowMat=as(matrix(estiResults$CIvecLow,nrow=nPredics),"sparseMatrix")
  results$CIUpMat=as(matrix(estiResults$CIvecUp,nrow=nPredics),"sparseMatrix")
  rm(estiResults)

  estByCovList=list()
  count=0
  for(i in 1:nTestCov){
    sigTaxaPosition=which(results$selecMatAll[i,]!=0)
    nrow=length(sigTaxaPosition)
    if(nrow==0)next
    count=count+1
    ncol=3
    estByCovMat=matrix(NA,nrow=nrow,ncol=ncol)

    for(j in 1:nrow){
      for(k in 1:ncol){
        if(k==1)estByCovMat[j,k]=results$betaMat[i,sigTaxaPosition[j]]
        if(k==2)estByCovMat[j,k]=results$CILowMat[i,sigTaxaPosition[j]]
        if(k==3)estByCovMat[j,k]=results$CIUpMat[i,sigTaxaPosition[j]]
      }
    }

    rownames(estByCovMat)=microbName[results$selecMatAll[i,]!=0]
    colnames(estByCovMat)=c("Beta","LowBound95%CI","UpBound95%CI")
    estByCovList[[count]]=estByCovMat
    names(estByCovList)[i]=testCovInOrder[i]
    rm(estByCovMat)
  }

  results$estByCovList=estByCovList
  rm(estByCovList)

  SigCovByTaxaList=list()
  count=0
  for(i in 1:nTaxa){
    sigCov=which(results$selecMatAll[,i]!=0)
    if(length(sigCov)==0)next
    count=count+1
    SigCovByTaxaList[[count]]=testCovInOrder[sigCov]
    names(SigCovByTaxaList)[[count]]=microbName[i]
  }
  results$SigCovByTaxaList=SigCovByTaxaList
  rm(SigCovByTaxaList,microbName)

  results$reguMethod=reguMethod
  results$independence=independence
  results$nSub=nSub
  results$nTaxa=nTaxa
  results$nPredics=nPredics
  return(results)

  rm(data)

  # return results
  results$posFoldRange=posFoldRange
  results$negFoldRange=negFoldRange
  results$taxaCor=taxaCor
  results$truSNR=truSNR
  results$predCor=predCor
  results$ConsBaseYesNo=ConsBaseYesNo
  results$scenario=scenario
  results$logNormMean=logNormMean
  results$logNormSD=logNormSD
  results$coefSparsity=coefSparsity
  results$baseDataSparsity=baseDataSparsity
  results$reguMethod=reguMethod
  results$nSub=nSub
  results$nTaxaTru=nTaxaTru
  results$nTaxa=nTaxa
  results$nPredics=nPredics
  results$nPredicsWithInt=nPredicsWithInt
  results$nPredTru=nPredTru
  results$nAssociations=nAssociations
  results$independence=independence
  results$identity=identity
  results$echoPropDiff=echoPropDiff
  results$seqDepthType=seqDepthType
  results$seqDepthDeno=seqDepthDeno
  results$echoAbundanDist=echoAbundanDist
  results$echoGammaPar=echoGammaPar
  results$echoEffectsDist=echoEffectsDist
  results$echoEffects=echoEffects
  results$fwerRate=fwerRate
  results$doPermut=doPermut
  results$nPermu=nPermu
  return(results)
}

