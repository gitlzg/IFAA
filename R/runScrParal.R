
##' @export


runScrParal=function(
  method=c("mcp"),
  data,
  testCovInd,
  testCovInOrder,
  testCovInNewNam,
  nRef,
  paraJobs,
  doPermut,
  permutY=F,
  nPermu,
  refTaxa=NULL,
  allFunc=allFunc,
  refReadsThresh,
  SDThresh,
  SDquantilThresh,
  balanceCut,
  Mprefix,
  covsPrefix,
  binPredInd,
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

  results$goodRefTaxaCandi=basicInfo$goodRefTaxaCandi

  nSub=basicInfo$nSub
  nSubNoReads=basicInfo$nSubNoReads
  maxTaxaNameNum=basicInfo$maxTaxaNameNum
  rm(basicInfo)
  gc()

  nNorm=nTaxa-1
  nAlphaNoInt=nPredics*nNorm
  nAlphaSelec=nPredics*nTaxa

  # make reference taxa list
  if(length(refTaxa)==0){
    if(nRef>length(taxaNames)){
      stop("Error: number of random reference taxa can not be larger than the number
           of candidate reference taxa. Try lower nRef")
    }
    if(length(seed)>0){
      set.seed(as.numeric(seed))
    }
    refTaxa=sample(taxaNames,nRef)
    }

  # overwrite nRef if the reference taxon is specified
  nRef=length(refTaxa)
  #
  ## run original data screen
  #
  screen1=originDataScreen(data=data,testCovInd=testCovInd,
                           nRef=nRef,refTaxa=refTaxa,paraJobs=paraJobs,
                           method=method,allFunc=allFunc,Mprefix=Mprefix,
                           covsPrefix=covsPrefix,
                           binPredInd=binPredInd,seed=seed)

  results$countOfSelecForAPred=screen1$countOfSelecForAPred
  yTildLongList=screen1$yTildLongList
  results$testCovCountMat=screen1$testCovCountMat
  rm(screen1)
  gc()

  #
  ## start to run permutation screen
  #
  if(doPermut){
    startT=proc.time()[3]
    cat("start to run permutation","\n")
    # permut the exposure variable
    permutOrder=lapply(rep(nSub,nPermu),sample)

    screenStartTime = proc.time()[3]
    refResu0=as(matrix(0,nrow=nPredics*nTaxa,ncol=1),"sparseMatrix")

    EName=testCovInNewNam
    EVar=data[,EName,drop=F]

    totNumOfLoops=nRef*nPermu

    if(length(paraJobs)==0){
      availCores=availableCores()
      if(is.numeric(availCores))paraJobs=max(1,availableCores()-2)
      if(!is.numeric(availCores))paraJobs=1
    }
    c2 <- snow::makeCluster(paraJobs)

    snow::clusterExport(c2, allFunc)

    if(length(seed)>0){
      snow::clusterSetupRNGstream(cl=c2,seed=as.numeric(seed)+10^2)
    }
    doSNOW::registerDoSNOW(c2)

    refResu=foreach (i=1:totNumOfLoops,.multicombine=T,
                     .packages=c("picasso","glmnet","expm","doSNOW","snow","foreach","Matrix"),
                     .errorhandling="pass") %dopar% {
                       #for(j in 1:nRef){
                       permut.i=1+(i-1)%%nPermu
                       permutX1=EVar[permutOrder[[permut.i]],,drop=F]
                       newData=cbind(data[,!(colnames(data)%in%EName)],permutX1)
                       rm(permutX1)
                       colnames(newData)[(ncol(newData)-length(EName)+1):ncol(newData)]=EName

                       ref.i=1+floor((i-1)/nPermu)
                       ii=which(taxaNames==refTaxa[ref.i])

                       xLong.i=dataRecovTrans(data=newData,ref=refTaxa[ref.i],
                                              Mprefix=Mprefix,covsPrefix=covsPrefix,xOnly=T)
                       xLongTild.i=xLong.i$xTildalong
                       rm(newData,xLong.i)
                       gc()
                       if(method=="lasso") {
                         Penal.i=runGlmnet(x=xLongTild.i,y=yTildLongList[[ref.i]],nPredics=nPredics)
                       }
                       if(method=="mcp") {
                         Penal.i=runPicasso(x=xLongTild.i,y=yTildLongList[[ref.i]],nPredics=nPredics,
                                            method="mcp",permutY=permutY)
                       }
                       rm(xLongTild.i)
                       BetaNoInt.i=as(Penal.i$betaNoInt,"sparseVector")
                       rm(Penal.i)
                       gc()
                       selection.i=as(rep(0,nAlphaSelec),"sparseVector")
                       if (ii==1){
                         selection.i[-seq(1,nPredics)]=as(BetaNoInt.i!=0,"sparseVector")
                       }
                       if (ii==nTaxa) {
                         selection.i[-seq((nAlphaSelec-nPredics+1),nAlphaSelec)]=as(BetaNoInt.i!=0,"sparseVector")
                       }
                       if ((ii>1) & (ii<nTaxa)) {
                         selection.i[1:(nPredics*(ii-1))]=as(BetaNoInt.i[1:(nPredics*(ii-1))]!=0,"sparseVector")
                         selection.i[(nPredics*ii+1):nAlphaSelec]=as(BetaNoInt.i[(nPredics*(ii-1)+1):nAlphaNoInt]!=0,"sparseVector")
                       }
                       rm(BetaNoInt.i)
                       # create return vector
                       recturnVec=as(rep(0,nAlphaSelec),"sparseVector")
                       recturnVec=selection.i
                       rm(selection.i)
                       return(recturnVec)
                     }
    rm(yTildLongList)
    snow::stopCluster(c2)
    gc()

    refResu<- lapply(refResu, as, "sparseMatrix")
    refResu=do.call(cbind, refResu)

    rm(data,permutOrder,EVar)

    endT=proc.time()[3]

    cat("Permutation done and took", (endT-startT)/60,"minutes","\n")

    # obtain the maximum vector
    permuColInd=1+seq(totNumOfLoops-1)%%nPermu
    permutResuMat=matrix(NA,nrow=nTaxa,ncol=nPermu)
    permutTestCovList=list()
    nTestCov=length(testCovInd)

    for(i in 1:nPermu){
      allCovCountMat.i=as(matrix(Matrix::rowSums(refResu[,(permuColInd==i),drop=F]),nrow=nPredics),"sparseMatrix")
      testCovCountMat.i=allCovCountMat.i[testCovInd,,drop=F]
      permutResuMat[,i]=Matrix::colSums(testCovCountMat.i)
      if(nTestCov>1)permutTestCovList[[i]]=testCovCountMat.i
    }
    rm(allCovCountMat.i,testCovCountMat.i)

    nTestCov=length(testCovInd)

    maxVec=rep(NA,nPermu)
    for(i in 1:nPermu){
      maxVec[i]=max(permutResuMat[,i])
    }

    results$nTestCov=nTestCov

    if(nTestCov>1){
      MaxMatTestCovByPermu=matrix(NA,nrow=nTestCov,ncol=nPermu)
      for(k in 1:nTestCov){
        for(i in 1:nPermu){
          MaxMatTestCovByPermu[k,i]=max(permutTestCovList[[i]][k,])
        }
      }
      rm(permutTestCovList)

      results$MaxMatTestCovByPermu=MaxMatTestCovByPermu
      rm(MaxMatTestCovByPermu)
    }

    # obtain the null binomial distribution for each taxa
    totSeleCount=Matrix::rowSums(refResu)
    rm(refResu)
    binomPar=totSeleCount/totNumOfLoops
    rm(totSeleCount,totNumOfLoops)

    results$maxVec=maxVec
    rm(maxVec)

    results$binomPar=binomPar
    rm(binomPar)
    results$nTaxa=nTaxa
    results$nPredics=nPredics
  }

  results$taxaNames=taxaNames
  rm(taxaNames)
  return(results)
  }
