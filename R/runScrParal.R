
#---------------------------------------------------------------------------------------
## parallel screen function to select independent reference taxa
#---------------------------------------------------------------------------------------

runScrParal=function(
  method=c("mcp"),
  data,
  testCovInd,
  testCovInOrder,
  testCovInNewNam,
  nRef,
  paraJobs,
  x1permut,
  nPermu,
  refTaxa,
  maxDimensionScr=0.8*434*10*10^4,
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
  nSub=basicInfo$nSub
  predNames=basicInfo$predNames
  
  results$goodRefTaxaCandi=basicInfo$goodRefTaxaCandi
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
  gc(reset=TRUE)
  
  screen1=originDataScreen(data=data,testCovInd=testCovInd,
                           nRef=nRef,refTaxa=refTaxa,
                           paraJobs=paraJobs,
                           method=method,allFunc=allFunc,Mprefix=Mprefix,
                           covsPrefix=covsPrefix,
                           binPredInd=binPredInd,
                           standardize=standardize,
                           sequentialRun=sequentialRun,
                           seed=seed)
  
  maxMemUsedInMb=sum(gc()[,6])
  message("maximum memory used after phase 1a: ",maxMemUsedInMb," Mb")
  
  results$countOfSelecForAPred=screen1$countOfSelecForAPred
  yTildLongList=screen1$yTildLongList
  results$testCovCountMat=screen1$testCovCountMat
  Alist=screen1$Alist
  Wlist=screen1$Wlist
  twoList=screen1$twoList
  lengthTwoList=screen1$lengthTwoList
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
    startT=proc.time()[3]
    message("start to run permutation")
    
    # permut the exposure variable
    if(length(seed)>0)set.seed(as.numeric(seed)+10^6)
    
    permutOrder=lapply(rep(nSub,nPermu),sample)
    if(!x1permut){
      residuPermuOrder=lapply(rep(length(residuList[[1]]),nPermu),sample)
    }
    
    screenStartTime = proc.time()[3]
    
    EName=testCovInNewNam
    EVar=data[,EName,drop=FALSE]
    
    totNumOfLoops=nRef*nPermu
    
    if(length(paraJobs)==0){
      availCores=availableCores()
      if(is.numeric(availCores))paraJobs=max(1,availableCores()-2)
      if(!is.numeric(availCores))paraJobs=1
    }
    
    if(!sequentialRun){
      message(paraJobs, " parallel jobs are registered for the permutation analysis in Phase 1b")
    }
    
    batch=paraJobs
    forLoopN=ceiling(totNumOfLoops/batch)

    gc(reset=TRUE)
    
    for (jj in 1:forLoopN){
      cl <- parallel::makeCluster(paraJobs)

      parallel::clusterExport(cl=cl, varlist=allFunc,envir=parent.env(environment()))
      doParallel::registerDoParallel(cl)
      if(length(seed)>0){
        set.seed(as.numeric(seed)+jj+10^6)
        parallel::clusterSetRNGStream(cl=cl,(as.numeric(seed)+jj+10^7))
      }
      if(sequentialRun){foreach::registerDoSEQ()}
      
      if(forLoopN>1 & jj<forLoopN){
        refResu.j=foreach(i=((jj-1)*batch+1):(jj*batch),
                          .multicombine=T,
                          .packages=c("picasso","Matrix"),
                          .errorhandling="pass") %dopar% {
            # for(i in ((jj-1)*paraJobs+1):(jj*paraJobs)){

                            permut.i=1+(i-1)%%nPermu
                            ref.i=1+floor((i-1)/nPermu)
                            ii=which(taxaNames==refTaxa[ref.i])

                            if(x1permut){
                              permutX1=EVar[permutOrder[[permut.i]],,drop=FALSE]
                              newData=data
                              newData[,EName]=permutX1
                              rm(permutX1)
                              
                              xLong.i=dataRecovTrans(data=newData,ref=refTaxa[ref.i],
                                                     Mprefix=Mprefix,covsPrefix=covsPrefix,xOnly=TRUE)
                              xLongTild.i=xLong.i$xTildalong
                              rm(newData,xLong.i)
                              
                              gc()
                              
                              maxSubSamplSiz=floor(maxDimensionScr/ncol(xLongTild.i))
                              nToSamplFrom=nrow(xLongTild.i)
                              
                              if(method=="mcp") {
                                subSamplK=ceiling(nToSamplFrom/maxSubSamplSiz)
                                if(subSamplK==1)maxSubSamplSiz=nToSamplFrom
                                
                                nRuns=ceiling(subSamplK/3)
                           
                                for (k in 1:nRuns){
                                  rowToKeep=sample(nToSamplFrom,maxSubSamplSiz)
                                  x=as((xLongTild.i[rowToKeep,]),"sparseMatrix")
                                  y=as((yTildLongList[[ref.i]][rowToKeep]),"sparseVector")
                                  Penal.i=runPicasso(x=x,y=y,
                                                     nPredics=nPredics,
                                                     method="mcp",permutY=FALSE,
                                                     standardize=standardize,
                                                     seed=seed,seedi=i)
                                  
                                  BetaNoInt.k=as((0+(Penal.i$betaNoInt!=0)),"sparseVector")
                                  rm(Penal.i)
                                  gc()
                                  if(k==1)BetaNoInt.i=BetaNoInt.k
                                  if(k>1)BetaNoInt.i=BetaNoInt.i+BetaNoInt.k
                                  rm(BetaNoInt.k)
                                }
                                rm(k,x,y,xLongTild.i)
                                gc()
                                BetaNoInt.i=BetaNoInt.i/nRuns
                              }
                              gc()
                            }else{
                              resid=residuList[[ref.i]][residuPermuOrder[[permut.i]]]
                              yTildLong.i=xBetaList[[ref.i]]+resid
                              rm(resid)

                              maxSubSamplSiz=floor(maxDimensionScr/ncol(xTildLong))
                              nToSamplFrom=nrow(xTildLong)
                              if(method=="mcp") {
                                subSamplK=ceiling(nToSamplFrom/maxSubSamplSiz)
                                if(subSamplK==1)maxSubSamplSiz=nToSamplFrom

                                nRuns=ceiling(subSamplK/3)

                                for (k in 1:nRuns){
                                  rowToKeep=sample(nToSamplFrom,maxSubSamplSiz)
                                  x=as((xTildLong[rowToKeep,]),"sparseMatrix")
                                  y=as((yTildLong.i[rowToKeep]),"sparseVector")
                                  Penal.i=runPicasso(x=x,y=y,
                                                     nPredics=nPredics,
                                                     method="mcp",permutY=FALSE,
                                                     standardize=standardize,
                                                     seed=seed,seedi=i)

                                  BetaNoInt.k=as((0+(Penal.i$betaNoInt!=0)),"sparseVector")
                                  rm(Penal.i)
                                  gc()
                                  if(k==1)BetaNoInt.i=BetaNoInt.k
                                  if(k>1)BetaNoInt.i=BetaNoInt.i+BetaNoInt.k
                                  rm(BetaNoInt.k)
                                }
                                rm(k,x,y,xTildLong)
                                gc()
                                BetaNoInt.i=BetaNoInt.i/nRuns
                              }
                              rm(yTildLong.i)
                            }
                            selection.i=as(rep(0,nAlphaSelec),"sparseVector")
                            if (ii==1){
                              selection.i[-seq(1,nPredics)]=BetaNoInt.i
                            }
                            if (ii==nTaxa) {
                              selection.i[-seq((nAlphaSelec-nPredics+1),nAlphaSelec)]=BetaNoInt.i
                            }
                            if ((ii>1) & (ii<nTaxa)) {
                              selection.i[1:(nPredics*(ii-1))]=BetaNoInt.i[1:(nPredics*(ii-1))]
                              selection.i[(nPredics*ii+1):nAlphaSelec]=BetaNoInt.i[(nPredics*(ii-1)+1):nAlphaNoInt]
                            }
                            rm(BetaNoInt.i)
                            # create return vector
                            recturnVec=as(rep(0,nAlphaSelec),"sparseVector")
                            recturnVec=selection.i
                            rm(selection.i)
                            return(recturnVec)
                          }
        parallel::stopCluster(cl)
        gc()
        if(jj==1)refResu=refResu.j
        if(jj>1)refResu=do.call(c,list(refResu,refResu.j))
      }
      if(jj==forLoopN){
        refResu.j=foreach(i=((forLoopN-1)*batch+1):totNumOfLoops,
                          .multicombine=T,
                          .packages=c("picasso","Matrix"),
                          .errorhandling="pass") %dopar% {

        # for(i in ((forLoopN-1)*batch+1):totNumOfLoops){
                            
                            permut.i=1+(i-1)%%nPermu
                            ref.i=1+floor((i-1)/nPermu)
                            ii=which(taxaNames==refTaxa[ref.i])
                            
                            if(x1permut){
                              permutX1=EVar[permutOrder[[permut.i]],,drop=FALSE]
                              newData=data
                              newData[,EName]=permutX1
                              rm(permutX1)
                              
                              xLong.i=dataRecovTrans(data=newData,ref=refTaxa[ref.i],
                                                     Mprefix=Mprefix,covsPrefix=covsPrefix,xOnly=TRUE)
                              xLongTild.i=xLong.i$xTildalong
                              rm(newData,xLong.i)

                              gc()

                              maxSubSamplSiz=floor(maxDimensionScr/ncol(xLongTild.i))
                              nToSamplFrom=nrow(xLongTild.i)
                              if(method=="mcp") {
                                subSamplK=ceiling(nToSamplFrom/maxSubSamplSiz)
                                if(subSamplK==1)maxSubSamplSiz=nToSamplFrom
                                nRuns=ceiling(subSamplK/3)

                                for (k in 1:nRuns){
                                  rowToKeep=sample(nToSamplFrom,maxSubSamplSiz)
                                  x=as((xLongTild.i[rowToKeep,]),"sparseMatrix")
                                  y=as((yTildLongList[[ref.i]][rowToKeep]),"sparseVector")
                                  Penal.i=runPicasso(x=x,y=y,
                                                     nPredics=nPredics,
                                                     method="mcp",permutY=FALSE,
                                                     standardize=standardize,
                                                     seed=seed,seedi=i)
                                  
                                  BetaNoInt.k=as((0+(Penal.i$betaNoInt!=0)),"sparseVector")
                                  rm(Penal.i)
                                  gc()
                                  if(k==1)BetaNoInt.i=BetaNoInt.k
                                  if(k>1)BetaNoInt.i=BetaNoInt.i+BetaNoInt.k
                                  rm(BetaNoInt.k)
                                }
                                rm(k,x,y,xLongTild.i)
                                gc()
                                BetaNoInt.i=BetaNoInt.i/nRuns
                              }
                              gc()
                            }else{
                              resid=residuList[[ref.i]][residuPermuOrder[[permut.i]]]
                              yTildLong.i=xBetaList[[ref.i]]+resid
                              rm(resid)
                              
                              maxSubSamplSiz=floor(maxDimensionScr/ncol(xTildLong))
                              nToSamplFrom=nrow(xTildLong)
                              if(method=="mcp") {
                                subSamplK=ceiling(nToSamplFrom/maxSubSamplSiz)
                                if(subSamplK==1)maxSubSamplSiz=nToSamplFrom

                                nRuns=ceiling(subSamplK/3)

                                for (k in 1:nRuns){
                                  rowToKeep=sample(nToSamplFrom,maxSubSamplSiz)
                                  x=xTildLong[rowToKeep,]
                                  y=yTildLong.i[rowToKeep]
                                  Penal.i=runPicasso(x=x,y=y,
                                                     nPredics=nPredics,
                                                     method="mcp",permutY=FALSE,
                                                     standardize=standardize,
                                                     seed=seed,seedi=i)

                                  BetaNoInt.k=as((0+(Penal.i$betaNoInt!=0)),"sparseVector")
                                  rm(Penal.i)
                                  gc()
                                  if(k==1)BetaNoInt.i=BetaNoInt.k
                                  if(k>1)BetaNoInt.i=BetaNoInt.i+BetaNoInt.k
                                  rm(BetaNoInt.k)
                                }
                                rm(k,x,y)
                                gc()
                                BetaNoInt.i=BetaNoInt.i/nRuns
                              }
                              rm(xTildLong,yTildLong.i)
                            }
                            
                            selection.i=as(rep(0,nAlphaSelec),"sparseVector")
                            if (ii==1){
                              selection.i[-seq(1,nPredics)]=BetaNoInt.i
                            }
                            if (ii==nTaxa) {
                              selection.i[-seq((nAlphaSelec-nPredics+1),nAlphaSelec)]=BetaNoInt.i
                            }
                            if ((ii>1) & (ii<nTaxa)) {
                              selection.i[1:(nPredics*(ii-1))]=BetaNoInt.i[1:(nPredics*(ii-1))]
                              selection.i[(nPredics*ii+1):nAlphaSelec]=BetaNoInt.i[(nPredics*(ii-1)+1):nAlphaNoInt]
                            }
                            rm(BetaNoInt.i)
                            # create return vector
                            recturnVec=as(rep(0,nAlphaSelec),"sparseVector")
                            recturnVec=selection.i
                            rm(selection.i)
                            return(recturnVec)
                          }
        parallel::stopCluster(cl)
        gc()
        if(forLoopN==1)refResu=refResu.j
        if(forLoopN>1)refResu=do.call(c,list(refResu,refResu.j))
      }
      
      if(jj>0 & (jj%%(ceiling(forLoopN/10))==0)){
        message(round(100*jj/forLoopN,0), " percent of phase 1b analysis has been done.")
      }
    }
    rm(Alist,Wlist,yTildLongList,twoList,lengthTwoList)
    gc()
    maxMemUsedInMb=sum(gc()[,6])
    message("maximum memory used after permutation: ",maxMemUsedInMb," Mb")

    refResu<- lapply(refResu, as, "sparseMatrix")
    refResu=do.call(cbind, refResu)
    
    rm(data,permutOrder,EVar)
    
    endT=proc.time()[3]
    
    message("Permutation analysis done and took ", round((endT-startT)/60,3)," minutes")
    
    # obtain the maximum vector
    permuColInd=1+seq(0,totNumOfLoops-1)%%nPermu
    permutResuMat=matrix(NA,nrow=nTaxa,ncol=nPermu)
    permutTestCovList=list()
    nTestCov=length(testCovInd)
    
    for(i in 1:nPermu){
      matrix.i.permu=refResu[,(permuColInd==i),drop=FALSE]
      vec.i=as(rep(0,nTaxa),"sparseVector")
      
      for(j in 1:nRef){
        vector.j.ref=matrix.i.permu[,j]
        matrix.i.j=as(matrix(vector.j.ref,nrow=nPredics),"sparseMatrix")
        testCovVec.j.ref=as((Matrix::colSums(matrix.i.j[testCovInd,,drop=FALSE])),"sparseVector")
        vec.i=vec.i+testCovVec.j.ref
      }
      
      permutResuMat[,i]=as.vector(vec.i)
    }
    rm(vec.i,matrix.i.permu,vector.j.ref,matrix.i.j,testCovVec.j.ref)
    
    if(nTestCov>1){
      for(i in 1:nPermu){
        matrix.i.permu=refResu[,(permuColInd==i),drop=FALSE]
        allCovCountMat.i=as(matrix(Matrix::rowSums(matrix.i.permu),nrow=nPredics),"sparseMatrix")
        permutTestCovList[[i]]=allCovCountMat.i[testCovInd,,drop=FALSE]
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

  results$taxaNames=taxaNames
  rm(taxaNames)
  return(results)
}

