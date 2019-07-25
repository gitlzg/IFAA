
##' @export


runScrParal=function(
  method=c("mcp"),
  data,
  testCovInd,
  nRef,
  paraJobs,
  doPermut,
  permutY=F,
  nPermu,
  refTaxa=NULL,
  independence,
  identity,
  allFunc=allFunc,
  refReadsThresh,
  SDThresh,
  SDquantilThresh,
  balanceCut,
  Mprefix,
  covsPrefix,
  binPredInd){

  results=list()

  # load data info
  basicInfo=dataInfo(data=data,refReadsThresh=refReadsThresh,
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
  EName=paste0(covsPrefix,1)
  EVar=data[,EName]
  EVarPermut=sample(EVar)
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
    set.seed(3)
    refTaxa=sample(taxaNames,nRef)
    }

  # overwrite nRef if the reference taxon is specified
  nRef=length(refTaxa)
  #
  ## run original data screen
  #
  screen1=originDataScreen(data=data,testCovInd=testCovInd,
                           nRef=nRef,refTaxa=refTaxa,paraJobs=paraJobs,
                           independence=independence,identity=identity,
                           method=method,allFunc=allFunc,Mprefix=Mprefix,
                           covsPrefix=covsPrefix,
                           binPredInd=binPredInd)

  results$countOfSelecForAPred=screen1$countOfSelecForAPred
  results$testCovCountMat=screen1$testCovCountMat
  results$MCPoriginSuccCV=screen1$MCPoriginSuccCV
  results$nCV=screen1$nCV
  rm(screen1)

  gc()

  #
  ## start to run permutation screen
  #
  if(doPermut){
    cat("start to run permutation","\n")
    # permut the exposure variable
    set.seed(4)
    permutOrder=lapply(rep(nSub,nPermu),sample)

    totNumOfLoops=nPermu*nRef
    screenStartTime = proc.time()[3]
    if(length(paraJobs)==0){
      availCores=availableCores()
      if(is.numeric(availCores))paraJobs=max(1,availableCores()-1)
      if(!is.numeric(availCores))paraJobs=1
    }

    c2 <- makeCluster(paraJobs)

    clusterExport(c2, varlist=allFunc)

    clusterSetupRNGstream(cl=c2,seed=2)

    registerDoSNOW(c2)

    refResu=foreach (i=1:totNumOfLoops,.multicombine=T,
                     .packages=c("picasso","glmnet","expm","doSNOW","snow","foreach","Matrix"),
                     .errorhandling="pass") %dopar% {
                       permut.i=1+i%%nPermu
                       permutX1=EVar[permutOrder[[permut.i]]]
                       newData=cbind(data[,!(colnames(data)%in%EName)],permutX1)
                       rm(permutX1)
                       colnames(newData)[ncol(newData)]=EName
                       ref.i=1+floor((i-1)/nPermu)
                       ii=which(taxaNames==refTaxa[ref.i])
                       if(permutY){
                         dataForEst=dataRecovTrans(data=data,independence=independence,
                                                   identity=identity,ref=refTaxa[ref.i],
                                                   Mprefix=Mprefix,covsPrefix=covsPrefix,
                                                   binPredInd=binPredInd)
                       }
                       else {
                         dataForEst=dataRecovTrans(data=newData,independence=independence,
                                                   identity=identity,ref=refTaxa[ref.i],
                                                   Mprefix=Mprefix,covsPrefix=covsPrefix,
                                                   binPredInd=binPredInd)
                       }
                       x=as(dataForEst$xTildalong,"sparseMatrix")
                       y=as(dataForEst$UtildaLong,"sparseVector")
                       rm(dataForEst)
                       gc()
                       if(method=="lasso") {
                         Penal.i=runGlmnet(x=x,y=y,nPredics=nPredics)
                         MCPpermuSuccCV=0
                       }

                       if(method=="mcp") {
                         Penal.i=runPicasso(x=x,y=y,nPredics=nPredics,method="mcp",
                                            permutY=permutY,allFunc=allFunc)
                         MCPpermuSuccCV=Penal.i$successCV
                       }
                       rm(x,y)
                       gc()

                       if(is.null(Penal.i)) {
                         Penal.i$betaNoInt=as(rep(0,nNorm),"sparseVector")
                       }
                       BetaNoInt.i=as(Penal.i$betaNoInt,"sparseVector")
                       rm(Penal.i)
                       gc()

                       selection.i=as(rep(0,nAlphaSelec),"sparseVector")

                       if (ii==1){
                         selection.i[-seq(1,nPredics)]=as(BetaNoInt.i!=0,"sparseVector")
                       }
                       if (ii==nTaxa) {
                         selection.i[-seq((nAlphaSelec-nPredics-1),nAlphaSelec)]=as(BetaNoInt.i!=0,"sparseVector")
                       }
                       if ((ii>1) & (ii<nTaxa)) {
                         selection.i[1:(nPredics*(ii-1))]=as(BetaNoInt.i[1:(nPredics*(ii-1))]!=0,"sparseVector")
                         selection.i[(nPredics*ii+1):nAlphaSelec]=as(BetaNoInt.i[(nPredics*(ii-1)+1):nAlphaNoInt]!=0,"sparseVector")
                       }
                       rm(BetaNoInt.i)
                       gc()

                       # create return vector
                       recturnVec=as(rep(0,(nAlphaSelec+1)),"sparseVector")
                       recturnVec[1]=MCPpermuSuccCV
                       recturnVec[-1]=selection.i
                       rm(selection.i)
                       return(recturnVec)
                     }
    registerDoSEQ()
    stopCluster(c2)
    rm(data,permutOrder,EVar)
    gc()

    refResu<- lapply(refResu, as, "sparseMatrix")
    refResu=do.call(cbind, refResu)

    MCPpermuSuccCV=t(matrix(refResu[1,],nrow=nPermu,ncol=nRef))

    # obtain the maximum vector
    permuColInd=1+seq(totNumOfLoops)%%nPermu
    permutResuMat=matrix(NA,nrow=nTaxa,ncol=nPermu)
    permutTestCovList=list()
    nTestCov=length(testCovInd)
    for(i in 1:nPermu){
      allCovCountMat=as(matrix(Matrix::rowSums(refResu[-1,(permuColInd==i),drop=F]),nrow=nPredics),"sparseMatrix")
      testCovCountMat=allCovCountMat[testCovInd,,drop=F]
      permutResuMat[,i]=Matrix::colSums(testCovCountMat)
      if(nTestCov>1)permutTestCovList[[i]]=testCovCountMat
    }
    rm(allCovCountMat,testCovCountMat)

    nTestCov=length(testCovInd)

    maxVec=rep(NA,nPermu)
    for(i in 1:nPermu){
      maxVec[i]=max(permutResuMat[,i])
    }

    results$nTestCov=nTestCov

    if(nTestCov>1){
      MatTestCovByTaxa=matrix(NA,nrow=nTestCov,ncol=nPermu)
      permutTestCovByTaxa=list()
      for(k in 1:nTaxa){
        for(i in 1:nPermu){
          MatTestCovByTaxa[,i]=permutTestCovList[[i]][,k]
        }
        permutTestCovByTaxa[[k]]=MatTestCovByTaxa
      }
      rm(MatTestCovByTaxa)

      results$permutTestCovByTaxa=permutTestCovByTaxa
      rm(permutTestCovByTaxa)
    }

    # obtain the null binomial distribution for each taxa
    totSeleCount=Matrix::rowSums(refResu[-1,,drop=F])
    rm(refResu)
    binomPar=totSeleCount/totNumOfLoops
    rm(totSeleCount,totNumOfLoops)

    results$maxVec=maxVec
    rm(maxVec)

    results$MCPpermuSuccCV=MCPpermuSuccCV
    rm(MCPpermuSuccCV)
    results$binomPar=binomPar
    rm(binomPar)
    results$nTaxa=nTaxa
    results$nPredics=nPredics
  }
  results$taxaNames=taxaNames
  rm(taxaNames)

  return(results)
}
