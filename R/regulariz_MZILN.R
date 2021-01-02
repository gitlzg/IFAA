

Regulariz_MZILN=function(
  data,
  testCovInd,
  testCovInOrder,
  microbName,
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

  newRefTaxNam=taxaNames[microbName%in%refTaxa]
  nRef=length(refTaxa)
  refTaxa_reOrder=microbName[microbName%in%refTaxa]

  reguResu=originDataScreen(method=reguMethod,data=data,testCovInd=testCovInd,
                            paraJobs=paraJobs,refTaxa=newRefTaxNam,standardize=standardize,
                            sequentialRun=sequentialRun,allFunc=allFunc,Mprefix=Mprefix,
                            covsPrefix=covsPrefix,binPredInd=binaryInd,seed=seed)

  results$SelecAllRefTaxaPred=reguResu$scr1ResuSelec
  rm(reguResu)

  message("Reference taxa are: ",refTaxa)

  betaMatList=list()
  CILowMatList=list()
  CIUpMatList=list()
  results$estiList=list()
  for(iii in 1:nRef){
    time11=proc.time()[3]
    originTaxNam=refTaxa_reOrder[iii]
    results$estiList[[originTaxNam]]=bootResuHDCI(data=data,
                             refTaxa=newRefTaxNam[iii],originRefTaxNam=originTaxNam,
                             bootB=bootB,bootLassoAlpha=bootLassoAlpha,
                             binPredInd=binaryInd,covsPrefix=covsPrefix,Mprefix=Mprefix,
                             paraJobs=paraJobs,
                             standardize=standardize,seed=seed)

    time12=proc.time()[3]
    message("Estimation done for the ", iii,"th reference taxon: ",refTaxa_reOrder[iii],
        " and it took ",round((time12-time11)/60,2)," minutes")

    estiResults=results$estiList[[refTaxa_reOrder[iii]]]

    betaMatList[[iii]]=as(matrix(estiResults$finalBetaEst.LPR,nrow=nPredics),"sparseMatrix")
    CILowMatList[[iii]]=as(matrix(estiResults$CIvecLow.LPR,nrow=nPredics),"sparseMatrix")
    CIUpMatList[[iii]]=as(matrix(estiResults$CIvecUp.LPR,nrow=nPredics),"sparseMatrix")
  }
  rm(estiResults)

  estByRefTaxaList=list()

  for(iii in 1:nRef){
    SelecRefTaxa.i=matrix(results$SelecAllRefTaxaPred[,iii],nrow=nPredics)
    estByCovList=list()

    for(i in 1:nTestCov){
      sigTaxaPosition=which(SelecRefTaxa.i[i,]!=0)
      nrow=length(sigTaxaPosition)
      if(nrow==0)next

      ncol=3
      estByCovMat=matrix(NA,nrow=nrow,ncol=ncol)

      for(j in 1:nrow){
        for(k in 1:ncol){
          if(k==1)estByCovMat[j,k]=betaMatList[[iii]][i,sigTaxaPosition[j]]
          if(k==2)estByCovMat[j,k]=CILowMatList[[iii]][i,sigTaxaPosition[j]]
          if(k==3)estByCovMat[j,k]=CIUpMatList[[iii]][i,sigTaxaPosition[j]]
        }
      }

      rownames(estByCovMat)=microbName[SelecRefTaxa.i[i,]!=0]
      colnames(estByCovMat)=c("Beta.LPR","LowB95%CI.LPR","UpB95%CI.LPR")

      estByCovList[[testCovInOrder[i]]]=estByCovMat
      rm(estByCovMat)
    }

    if(length(estByCovList)==0){
      results$estByRefTaxaList[[refTaxa_reOrder[iii]]][["estByCovList"]]="No significant assoication is identified."
    }else{
      results$estByRefTaxaList[[refTaxa_reOrder[iii]]][["estByCovList"]]=estByCovList
    }
  }


  rm(estByCovList)

  for(iii in 1:nRef){
    SelecRefTaxa.i=matrix(results$SelecAllRefTaxaPred[,iii],nrow=nPredics)
    SigCovByTaxaList=list()

    for(i in 1:nTaxa){
      sigCov=which(SelecRefTaxa.i[,i]!=0)
      if(length(sigCov)==0)next
      SigCovByTaxaList[[microbName[i]]]=testCovInOrder[sigCov]
    }
    if(length(SigCovByTaxaList)==0){
      results$estByRefTaxaList[[refTaxa_reOrder[iii]]][["SigCovByTaxaList"]]="No significant assoication is identified."
    }else{
      results$estByRefTaxaList[[refTaxa_reOrder[iii]]][["SigCovByTaxaList"]]=SigCovByTaxaList
    }
  }

  rm(SigCovByTaxaList,microbName)

  results$reguMethod=reguMethod
  results$nSub=nSub
  results$nTaxa=nTaxa
  results$nPredics=nPredics
  rm(data)
  # return results
  return(results)
}

