
#---------------------------------------------------------------------------------------
## Original screen function to select independent reference taxa
#---------------------------------------------------------------------------------------
originDataScreen=function(
  method,
  data,
  testCovInd,
  nRef,
  paraJobs,
  refTaxa,
  maxDimensionScr=0.8*434*10*10^4,
  #maxDimensionScr=10^5,
  standardize,
  sequentialRun,
  allFunc,
  Mprefix,
  covsPrefix,
  binPredInd,
  seed){
  
  results=list()
  
  # load data info
  basicInfo=dataInfo(data=data,Mprefix=Mprefix,
                     covsPrefix=covsPrefix,
                     binPredInd=binPredInd)
  
  taxaNames=basicInfo$taxaNames
  nTaxa=basicInfo$nTaxa
  nPredics=basicInfo$nPredics
  rm(basicInfo)
  gc()
  
  nNorm=nTaxa-1
  nAlphaNoInt=nPredics*nNorm
  nAlphaSelec=nPredics*nTaxa
  
  countOfSelec=rep(0,nAlphaSelec)
  resultsByRefTaxon=list()
  
  # overwrite nRef if the reference taxon is specified
  nRef=length(refTaxa)
  
  startT1=proc.time()[3]
  message("start phase 1a")
  if(length(paraJobs)==0){
    availCores=availableCores()
    if(is.numeric(availCores))paraJobs=max(1,availableCores()-2)
    if(!is.numeric(availCores))paraJobs=1
  }
  batch=paraJobs
  forLoopN=ceiling(nRef/batch)

  if(!sequentialRun){
    message(paraJobs, " parallel jobs are registered for analyzing ", nRef, " reference taxa in Phase 1a")
  }
  
  for (jj in 1:forLoopN){

    cl<-parallel::makeCluster(paraJobs)
    
    parallel::clusterExport(cl=cl, varlist=allFunc,envir=parent.env(environment()))
    doParallel::registerDoParallel(cl)
    
    if(length(seed)>0){
      set.seed(as.numeric(seed)+jj+10^4)
      parallel::clusterSetRNGStream(cl=cl,(as.numeric(seed)+jj+10^5))
    }
    
    if(sequentialRun){foreach::registerDoSEQ()}
    
    #start parallel computing
    if(forLoopN>1 & jj<forLoopN){
      scr1Resu.j=foreach(i=((jj-1)*batch+1):(jj*batch),
                         .multicombine=T,
                         .packages=c("picasso","Matrix"),
                         .errorhandling="pass") %dopar% {
                           
                           # for(i in 1:((jj-1)*paraJobs+1):(jj*paraJobs)){

                           ii=which(taxaNames==refTaxa[i])
                           dataForEst=dataRecovTrans(data=data,ref=refTaxa[i],Mprefix=Mprefix,
                                                     covsPrefix=covsPrefix)
                           
                           xTildLongTild.i=dataForEst$xTildalong
                           yTildLongTild.i=dataForEst$UtildaLong
                           rm(dataForEst)
                           gc()

                           maxSubSamplSiz=floor(maxDimensionScr/ncol(xTildLongTild.i))
                           nToSamplFrom=nrow(xTildLongTild.i)
                           
                           if(method=="mcp") {
                             subSamplK=ceiling(nToSamplFrom/maxSubSamplSiz)
                             if(subSamplK==1)maxSubSamplSiz=nToSamplFrom

                             nRuns=ceiling(subSamplK/3)

                             for (k in 1:nRuns){
                               rowToKeep=sample(nToSamplFrom,maxSubSamplSiz)
                               x=as((xTildLongTild.i[rowToKeep,]),"sparseMatrix")
                               y=as(yTildLongTild.i[rowToKeep],"sparseVector")
                               Penal.i=runPicasso(x=x,y=y,
                                                  nPredics=nPredics,
                                                  method="mcp",permutY=FALSE,
                                                  standardize=standardize,
                                                  seed=seed,seedi=i)
                               
                               BetaNoInt.k=as((0+(Penal.i$betaNoInt!=0)),"sparseVector")
                               rm(Penal.i)
                               if(k==1)BetaNoInt.i=BetaNoInt.k
                               if(k>1)BetaNoInt.i=BetaNoInt.i+BetaNoInt.k
                               rm(BetaNoInt.k)
                             }
                             rm(k,x,y,xTildLongTild.i)
                             gc()
                             BetaNoInt.i=BetaNoInt.i/nRuns
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
                           recturnlist=list()
                           recturnlist[[1]]=selection.i
                           recturnlist[[2]]=yTildLongTild.i
                           rm(selection.i,yTildLongTild.i)
                           return(recturnlist)
                         }
      parallel::stopCluster(cl)
      gc()
      if(jj==1)scr1Resu=scr1Resu.j
      if(jj>1)scr1Resu=do.call(c,list(scr1Resu,scr1Resu.j))
    }
    
    if(jj==forLoopN){
      scr1Resu.j=foreach(i=((forLoopN-1)*batch+1):nRef,
                         .multicombine=T,
                         .packages=c("picasso","Matrix"),
                         .errorhandling="pass") %dopar% {
                           
                           # for(i in ((forLoopN-1)*paraJobs+1):nRef){
                           
                           ii=which(taxaNames==refTaxa[i])
                           dataForEst=dataRecovTrans(data=data,ref=refTaxa[i],Mprefix=Mprefix,
                                                     covsPrefix=covsPrefix)
                           
                           xTildLongTild.i=dataForEst$xTildalong
                           yTildLongTild.i=dataForEst$UtildaLong
                           rm(dataForEst)
                           gc()

                           maxSubSamplSiz=floor(maxDimensionScr/ncol(xTildLongTild.i))
                           nToSamplFrom=nrow(xTildLongTild.i)

                           if(method=="mcp") {
                             
                             subSamplK=ceiling(nToSamplFrom/maxSubSamplSiz)
                             if(subSamplK==1)maxSubSamplSiz=nToSamplFrom

                             nRuns=ceiling(subSamplK/3)

                             for (k in 1:nRuns){
                               rowToKeep=sample(nToSamplFrom,maxSubSamplSiz)
                               x=as((xTildLongTild.i[rowToKeep,]),"sparseMatrix")
                               y=as((yTildLongTild.i[rowToKeep]),"sparseVector")
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
                             rm(k,x,y,xTildLongTild.i)
                             gc()
                             BetaNoInt.i=BetaNoInt.i/nRuns
                           }
                           
                           selection.i=as(rep(0,nAlphaSelec),"sparseVector")
                           if (ii==1){selection.i[-seq(1,nPredics)]=BetaNoInt.i
                           }
                           if (ii==nTaxa) {selection.i[-seq((nAlphaSelec-nPredics+1),nAlphaSelec)]=BetaNoInt.i
                           }
                           if ((ii>1) & (ii<nTaxa)) {
                             selection.i[1:(nPredics*(ii-1))]=BetaNoInt.i[1:(nPredics*(ii-1))]
                             selection.i[(nPredics*ii+1):nAlphaSelec]=BetaNoInt.i[(nPredics*(ii-1)+1):nAlphaNoInt]
                           }
                           rm(BetaNoInt.i)
                           # create return vector
                           recturnlist=list()
                           recturnlist[[1]]=selection.i
                           recturnlist[[2]]=yTildLongTild.i
                           rm(selection.i,yTildLongTild.i)
                           return(recturnlist)
                         }
      parallel::stopCluster(cl)
      if(forLoopN==1)scr1Resu=scr1Resu.j
      if(forLoopN>1)scr1Resu=do.call(c,list(scr1Resu,scr1Resu.j))
      gc()
    }

    if(jj>0 & (jj%%(ceiling(forLoopN/10))==0)){
      message(round(100*jj/forLoopN,0), " percent of phase 1a analysis has been done")
    }
  }
  rm(data)
  endT=proc.time()[3]
  
  message("Phase 1a done and took ",round((endT-startT1)/60,3)," minutes")
  
  selecList=list()
  for(i in 1:nRef){
    selecList[[i]]=scr1Resu[[i]][[1]]
  }
  
  results$yTildLongList=list()
  for(i in 1:nRef){
    results$yTildLongList[[i]]=scr1Resu[[i]][[2]]
  }

  rm(scr1Resu)

  selecList<- lapply(selecList, as, "sparseMatrix")
  scr1ResuSelec=do.call(cbind, selecList)
  rm(selecList)
  
  results$scr1ResuSelec=scr1ResuSelec
  
  # create count of selection for individual testCov
  countOfSelecForAllPred=as(matrix(Matrix::rowSums(scr1ResuSelec),nrow=nPredics),"sparseMatrix")
  testCovCountMat=countOfSelecForAllPred[testCovInd,,drop=FALSE]
  rm(scr1ResuSelec,testCovInd,countOfSelecForAllPred)

  # create overall count of selection for all testCov as a whole
  countOfSelecForAPred=as(matrix(Matrix::colSums(testCovCountMat),nrow=1),"sparseMatrix")
  gc()
  
  colnames(countOfSelecForAPred)=taxaNames
  rm(taxaNames)

  # return results
  results$testCovCountMat=testCovCountMat
  rm(testCovCountMat)
  results$countOfSelecForAPred=countOfSelecForAPred
  rm(countOfSelecForAPred)
  return(results)
}
