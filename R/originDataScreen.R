##' @export

originDataScreen=function(
  method,
  data,
  testCovInd,
  nRef,
  paraJobs,
  lambda=NULL,
  refTaxa,
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
  
  startT=proc.time()[3]
  cat("start Original screen","\n")
  if(length(paraJobs)==0){
    availCores=availableCores()
    if(is.numeric(availCores))paraJobs=max(1,availableCores()-2)
    if(!is.numeric(availCores))paraJobs=1
  }
  
  c1<-snow::makeCluster(paraJobs)
  
  if(!sequentialRun){
    cat(paraJobs, "parallel jobs are registered for analyzing", nRef, "reference taxa in Phase 1a.","\n")
  }
  
  snow::clusterExport(c1, allFunc)
  doSNOW::registerDoSNOW(c1)
  
  if(sequentialRun){foreach::registerDoSEQ()}
  
  startT1=proc.time()[3]
  cat("OriginDataScreen parallel setup took",startT1-startT,"seconds","\n")
  # start parallel computing
  scr1Resu=foreach(i=1:nRef,.multicombine=T,
                   .packages=c("picasso","expm","doSNOW","snow","foreach","Matrix"),
                   .errorhandling="pass") %dopar% {
                     
                     ii=which(taxaNames==refTaxa[i])
                     dataForEst=dataRecovTrans(data=data,ref=refTaxa[i],Mprefix=Mprefix,
                                               covsPrefix=covsPrefix)
                     xTildLongTild.i=dataForEst$xTildalong
                     yTildLongTild.i=dataForEst$UtildaLong
                     rm(dataForEst)
                     gc()
                     
                     if(method=="mcp") {
                       Penal.i=runPicasso(x=xTildLongTild.i,y=yTildLongTild.i,
                                          lambda=lambda,nPredics=nPredics,
                                          method="mcp",permutY=F,
                                          standardize=standardize,
                                          seed=seed,seedi=i)
                     }
                     rm(xTildLongTild.i)
                     
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
                     time4=proc.time()[3]
                     # create return vector
                     recturnlist=list()
                     recturnlist[[1]]=selection.i
                     recturnlist[[2]]=yTildLongTild.i
                     
                     rm(selection.i,yTildLongTild.i)
                     return(recturnlist)
                   }
  snow::stopCluster(c1)
  rm(data)
  gc()
  
  endT=proc.time()[3]
  
  cat("Original screen done and took",(endT-startT1)/60,"minutes","\n")
  
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
  testCovCountMat=countOfSelecForAllPred[testCovInd,,drop=F]
  rm(testCovInd,countOfSelecForAllPred)
  
  # create overall count of selection for all testCov as a whole
  countOfSelecForAPred=as(matrix(rep(0,nTaxa),nrow=1),"sparseMatrix")
  for (tax in 1:nTaxa){
    countMatForTaxni=scr1ResuSelec[(1+(tax-1)*nPredics):(tax*nPredics),,drop=F]
    totCountVecForTaxoni=Matrix::colSums(countMatForTaxni)
    countOfSelecForAPred[1,tax]=sum(totCountVecForTaxoni>0)
  }
  rm(tax,scr1ResuSelec,countMatForTaxni,totCountVecForTaxoni)
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


