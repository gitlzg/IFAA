
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
  x1permut,
  nPermu,
  refTaxa,
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
    refTaxa_extra=sample((taxaNames[!(taxaNames%in%refTaxa)]),(nRef-length(refTaxa)))
    refTaxa=c(refTaxa,refTaxa_extra)
    results$refTaxa=refTaxa
  }
  
  if(length(refTaxa)>=nRef){
    if(length(seed)>0){
      set.seed(as.numeric(seed))
    }
    refTaxa=sample(refTaxa,nRef)
    results$refTaxa=refTaxa
  }
  #
  ## run original data screen
  #
  screen1=originDataScreen(data=data,testCovInd=testCovInd,
                           nRef=nRef,refTaxa=refTaxa,paraJobs=paraJobs,
                           method=method,allFunc=allFunc,Mprefix=Mprefix,
                           covsPrefix=covsPrefix,
                           binPredInd=binPredInd,
                           standardize=standardize,
                           sequentialRun=sequentialRun,
                           seed=seed)
  
  results$countOfSelecForAPred=screen1$countOfSelecForAPred
  yTildLongList=screen1$yTildLongList
  results$testCovCountMat=screen1$testCovCountMat
  rm(screen1)
  gc()
  
  if(!x1permut){
    # generate Xbeta and residues for reduced model permutation
    XbetResi=XbetaAndResidu(data=data,testCovInd=testCovInd,
                            nRef=nRef,refTaxa=refTaxa,paraJobs=paraJobs,
                            method=method,allFunc=allFunc,Mprefix=Mprefix,
                            covsPrefix=covsPrefix,standardize=standardize,
                            binPredInd=binPredInd,
                            sequentialRun=sequentialRun,
                            seed=seed)
    xBetaList=XbetResi$xBetaList
    residuList=XbetResi$residuList
    xTildLong=XbetResi$xTildLong
    rm(XbetResi)
  }
  #
  ## start to run permutation screen
  #
  if(doPermut){
    startT=proc.time()[3]
    cat("start to run permutation","\n")
    
    # permut the exposure variable
    if(length(seed)>0)set.seed(as.numeric(seed)+10^6)
    
    permutOrder=lapply(rep(nSub,nPermu),sample)
    if(!x1permut){
      residuPermuOrder=lapply(rep(length(residuList[[1]]),nPermu),sample)
    }
    
    screenStartTime = proc.time()[3]
    
    EName=testCovInNewNam
    EVar=data[,EName,drop=F]
    
    totNumOfLoops=nRef*nPermu
    
    if(length(paraJobs)==0){
      availCores=availableCores()
      if(is.numeric(availCores))paraJobs=max(1,availableCores()-2)
      if(!is.numeric(availCores))paraJobs=1
    }
    
    c2 <- parallel::makeCluster(paraJobs)
    
    if(!sequentialRun){
      cat(paraJobs, "parallel jobs are registered for the permutation analysis in Phase 1b.","\n")
    }
    
    parallel::clusterExport(cl=c2,varlist=allFunc,envir=parent.env(environment()))
    doParallel::registerDoParallel(c2)
    
    if(sequentialRun){foreach::registerDoSEQ()}
    
    refResu=foreach (i=1:totNumOfLoops,.multicombine=T,
                     .packages=c("picasso","expm","foreach","Matrix"),
                     .errorhandling="pass") %dopar% {
                       #for(j in 1:nRef){
                       
                       permut.i=1+(i-1)%%nPermu
                       ref.i=1+floor((i-1)/nPermu)
                       ii=which(taxaNames==refTaxa[ref.i])
                       
                       if(x1permut){
                         permutX1=EVar[permutOrder[[permut.i]],,drop=F]
                         newData=data
                         newData[,EName]=permutX1
                         rm(permutX1)
                         
                         xLong.i=dataRecovTrans(data=newData,ref=refTaxa[ref.i],
                                                Mprefix=Mprefix,covsPrefix=covsPrefix,xOnly=T)
                         xLongTild.i=xLong.i$xTildalong
                         rm(newData,xLong.i)
                         gc()
                         
                         if(method=="mcp") {
                           Penal.i=runPicasso(x=xLongTild.i,y=yTildLongList[[ref.i]],nPredics=nPredics,
                                              method="mcp",permutY=permutY,
                                              standardize=standardize,
                                              seed=seed,seedi=i)
                         }
                         rm(xLongTild.i)
                       }else{
                         resid=residuList[[ref.i]][residuPermuOrder[[permut.i]]]
                         yTildLong.i=xBetaList[[ref.i]]+resid
                         rm(resid)
                         
                         if(method=="mcp") {
                           Penal.i=runPicasso(x=xTildLong,y=yTildLong.i,nPredics=nPredics,
                                              method="mcp",permutY=permutY,
                                              standardize=standardize,
                                              seed=seed,seedi=i)
                         }
                         rm(yTildLong.i)
                       }
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
    parallel::stopCluster(c2)
    gc()
    
    refResu<- lapply(refResu, as, "sparseMatrix")
    refResu=do.call(cbind, refResu)
    
    rm(data,permutOrder,EVar)
    
    endT=proc.time()[3]
    
    cat("Permutation done and took", (endT-startT)/60,"minutes","\n")
    
    # obtain the maximum vector
    permuColInd=1+seq(0,totNumOfLoops-1)%%nPermu
    permutResuMat=matrix(NA,nrow=nTaxa,ncol=nPermu)
    permutTestCovList=list()
    nTestCov=length(testCovInd)
    
    for(i in 1:nPermu){
      matrix.i.permu=refResu[,(permuColInd==i),drop=F]
      vec.i=as(rep(0,nTaxa),"sparseVector")
      
      for(j in 1:nRef){
        vector.j.ref=matrix.i.permu[,j]
        matrix.i.j=as(matrix(vector.j.ref,nrow=nPredics),"sparseMatrix")
        testCovVec.j.ref=as((Matrix::colSums(matrix.i.j[testCovInd,,drop=F])>0)+0,"sparseVector")
        vec.i=vec.i+testCovVec.j.ref
      }
      
      permutResuMat[,i]=as.vector(vec.i)
    }
    rm(vec.i,matrix.i.permu,vector.j.ref,matrix.i.j,testCovVec.j.ref)
    
    if(nTestCov>1){
      for(i in 1:nPermu){
        matrix.i.permu=refResu[,(permuColInd==i),drop=F]
        allCovCountMat.i=as(matrix(Matrix::rowSums(matrix.i.permu),nrow=nPredics),"sparseMatrix")
        permutTestCovList[[i]]=allCovCountMat.i[testCovInd,,drop=F]
      }
      rm(matrix.i.permu,allCovCountMat.i)
    }
    
    maxVec=rep(NA,nPermu)
    for(i in 1:nPermu){
      maxVec[i]=max(permutResuMat[,i])
    }
    rm(permutResuMat)
    
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

