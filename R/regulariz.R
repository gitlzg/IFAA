
##' @export



Regulariz=function(
  data,
  testCovInd,
  testCovInOrder,
  testCovInNewNam,
  microbName,
  nRef,
  nPermu,
  refTaxa,
  paraJobs,
  binaryInd,
  covsPrefix,
  Mprefix,
  reguMethod,
  fwerRate,
  boot,
  bootB,
  bootLassoAlpha,
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
  
  if(length(refTaxa)==0 & length(reguMethod)==1){
    # run MCP penalization
    
    regul.start.time = proc.time()[3]
    cat("Start Phase 1 association identification","\n")
    selectRegroup=getScrResu(data=data,testCovInd=testCovInd,
                             testCovInOrder=testCovInOrder,
                             testCovInNewNam=testCovInNewNam,nRef=nRef,
                             nPermu=nPermu,paraJobs=paraJobs,
                             refTaxa=NULL,fwerRate=fwerRate,
                             method=reguMethod,allFunc=allFunc,
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
    results$selecCountMatIndv=selectRegroup$selecCountMatIndv
    results$maxVec=selectRegroup$maxVec
    results$MaxMatTestCovByPermu=selectRegroup$MaxMatTestCovByPermu
    
    finalIndpRefTax=selectRegroup$finalIndpRefTax
    results$finalRefTaxonQualified=selectRegroup$refTaxonQualified
    
    results$goodIndpRefTaxLeastCount=selectRegroup$goodIndpRefTaxLeastCount
    results$goodIndpRefTaxWithCount=selectRegroup$goodIndpRefTaxWithCount
    results$goodIndpRefTaxFWERcut=selectRegroup$goodIndpRefTaxFWERcut
    
    results$goodIndpRefTaxFWERcutLeastCount=selectRegroup$goodIndpRefTaxFWERcutLeastCount
    results$goodIndpRefTaxFWERcut=selectRegroup$goodIndpRefTaxFWERcut
    results$goodRefTaxaCandi=selectRegroup$goodRefTaxaCandi
    
    results$taxaLessGoodCut=selectRegroup$taxaLessGoodCut
    results$taxaLessFWERCut=selectRegroup$taxaLessFWERCut
    results$refTaxa=selectRegroup$scrParal
    
    results$fwerCut=selectRegroup$fwerCut
    results$twoMeanUsed=selectRegroup$twoMeanUsed
    results$fwerCutIndv=selectRegroup$fwerCutIndv
    results$taxaNames=selectRegroup$taxaNames
    rm(selectRegroup)
    
    MCPExecuTime = (proc.time()[3] - regul.start.time)/60
    results$MCPExecuTime=MCPExecuTime
    cat("Phase 1 Associaiton identification is done and used", MCPExecuTime,"minutes","\n")
    
    results$finalizedBootRefTaxon=finalIndpRefTax
    cat("Final Reference Taxon is",results$finalizedBootRefTaxon,"\n")
    
    startT=proc.time()[3]
    cat("Start Phase 2 parameter estimation","\n")
    
    unestimableTaxa=c()
    qualifyData=data
    if(length(binaryInd)>0){
      qualifyData=data[rowSums(data[,results$taxaNames]>0)>=2,,drop=F]
      firstBinPredNam=paste0(covsPrefix,binaryInd)
      allBinPred=paste0(covsPrefix,binaryInd:nPredics)
      nBinPred=length(allBinPred)
      
      # find the pairs of binary preds and taxa for which the assocaiton is not identifiable
      AllTaxaNamesNoRefTax=results$taxaNames[!results$taxaNames%in%results$finalizedBootRefTaxon]
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
    TaxaNoReads=which(Matrix::colSums(qualifyData[,results$taxaNames])==0)
    rm(qualifyData)
    unestimableTaxa=unique(c(unestimableTaxa,results$taxaNames[TaxaNoReads]))
    results$unEstTaxaPos=which(results$taxaNames%in%unestimableTaxa)
    rm(TaxaNoReads,unestimableTaxa)
    
    estiResults=bootResuHDCI(data=data,refTaxa=results$finalizedBootRefTaxon,
                             bootB=bootB,bootLassoAlpha=bootLassoAlpha,
                             binPredInd=binaryInd,
                             covsPrefix=covsPrefix,Mprefix=Mprefix,
                             seed=seed)
    endT=proc.time()[3]
    cat("Phase 2 parameter estimation done and took",(endT-startT)/60,"minutes.","\n")
  }
  if(length(refTaxa)==1){
    CItime1 = proc.time()[3]
    refTaxaWithNewNam=paste0(Mprefix,which(microbName==refTaxa))
    results$finalizedBootRefTaxon=refTaxaWithNewNam
    
    cat("Final Reference Taxon is",results$finalRefTaxon,"\n")
    
    cat("Start association estimation","\n")
    estiResults=bootResuHDCI(data=data,refTaxa=refTaxaWithNewNam,
                             bootB=bootB,bootLassoAlpha=bootLassoAlpha,
                             binPredInd=binaryInd,
                             covsPrefix=covsPrefix,Mprefix=Mprefix,
                             seed=seed)
    
    cat("Association estimation done","\n")
    CItime2 = proc.time()[3]
    CIduration=(CItime2-CItime1)/60
    cat("Calculation of CI took", CIduration,"minutes","\n")
    
    results$betaMat=as(matrix(estiResults$finalBetaEst,nrow=nPredics),"sparseMatrix")
    results$CILowMat=as(matrix(estiResults$CIvecLow,nrow=nPredics),"sparseMatrix")
    results$CIUpMat=as(matrix(estiResults$CIvecUp,nrow=nPredics),"sparseMatrix")
    rm(estiResults)
    return(results)
  }
  
  results$betaMat=as(matrix(estiResults$finalBetaEst,nrow=nPredics),"sparseMatrix")
  results$CILowMat=as(matrix(estiResults$CIvecLow,nrow=nPredics),"sparseMatrix")
  results$CIUpMat=as(matrix(estiResults$CIvecUp,nrow=nPredics),"sparseMatrix")
  rm(estiResults)
  
  estByCovList=list()
  count=0
  
  for(i in 1:nTestCov){
    sigTaxaPosition=which(results$selecMatIndv[i,]!=0)
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
    
    rownames(estByCovMat)=microbName[results$selecMatIndv[i,]!=0]
    colnames(estByCovMat)=c("Beta","LowBound95%CI","UpBound95%CI")  
    
    estByCovList[[count]]=estByCovMat
    names(estByCovList)[count]=testCovInOrder[i]
    rm(estByCovMat)
  }
  
  if(length(estByCovList)==0){
    results$estByCovList="No significant assoication is identified. You may consider a higher FDR level."
  }else{
    results$estByCovList=estByCovList
  }
  rm(estByCovList)
  
  SigCovByTaxaList=list()
  count=0
  for(i in 1:nTaxa){
    sigCov=which(results$selecMatIndv[,i]!=0)
    if(length(sigCov)==0)next
    count=count+1
    SigCovByTaxaList[[count]]=testCovInOrder[sigCov]
    names(SigCovByTaxaList)[[count]]=microbName[i]
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
