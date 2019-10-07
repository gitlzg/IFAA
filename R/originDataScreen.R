##' @export


originDataScreen=function(
  method,
  data,
  testCovInd,
  nRef,
  paraJobs,
  lambda=NULL,
  refTaxa,
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
  
  snow::clusterExport(c1, allFunc)
  #cat("seed:",seed,"\n")
  
  #if(length(seed)>0){
  #snow::clusterSetupRNG (cl=c1, type = "RNGstream")
  #snow::clusterSetupRNGstream(cl=c1,seed=rep((as.numeric(seed)+10^5),6))
  #registerDoRNG((as.numeric(seed)+10^5))
  # }
  doSNOW::registerDoSNOW(c1)
  startT1=proc.time()[3]
  cat("OriginDataScreen parallel setup took",startT1-startT,"seconds","\n")
  # start parallel computing
  scr1Resu=foreach(i=1:nRef,.multicombine=T,
                   .packages=c("picasso","glmnet","expm","doSNOW","snow","foreach","Matrix"),
                   .errorhandling="pass") %dopar% {
                     #for(i in 1:nRef) {
                     
                     time1=proc.time()[3]
                     ii=which(taxaNames==refTaxa[i])
                     dataForEst=dataRecovTrans(data=data,ref=refTaxa[i],Mprefix=Mprefix,
                                               covsPrefix=covsPrefix)
                     xTildLongTild.i=dataForEst$xTildalong
                     yTildLongTild.i=dataForEst$UtildaLong
                     rm(dataForEst)
                     gc()
                     time2=proc.time()[3]
                     
                     preprosT=time2-time1
                     
                     if(method=="lasso") {
                       Penal.i=runGlmnet(x=xTildLongTild.ig,y=yTildLongTild.i,nPredics=nPredics)
                     }
                     if(method=="mcp") {
                       Penal.i=runPicasso(x=xTildLongTild.i,y=yTildLongTild.i,
                                          lambda=lambda,nPredics=nPredics,
                                          method="mcp",permutY=F,seed=seed,i=i)
                     }
                     rm(xTildLongTild.i)
                     
                     time3=proc.time()[3]
                     runpicasT=time3-time2
                     
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
                     #rm(BetaNoInt.i)
                     time4=proc.time()[3]
                     postProsT=time4-time3
                     # create return vector
                     recturnlist=list()
                     recturnlist[[1]]=selection.i
                     recturnlist[[2]]=yTildLongTild.i
                     recturnlist[[3]]=preprosT
                     recturnlist[[4]]=postProsT
                     recturnlist[[5]]=runpicasT
                     recturnlist[[6]]=BetaNoInt.i
                     rm(selection.i,yTildLongTild.i)
                     return(recturnlist)
                   }
  snow::stopCluster(c1)
  gc()
  
  endT=proc.time()[3]
  
  cat("Original screen done and took",(endT-startT1)/60,"minutes","\n")
  
  results$preprosT=vector()
  for(i in 1:nRef){
    results$preprosT[i]=scr1Resu[[i]][[3]]
  }
  
  results$postProsT=vector()
  for(i in 1:nRef){
    results$postProsT[i]=scr1Resu[[i]][[4]]
  }
  
  results$runpicasT=vector()
  for(i in 1:nRef){
    results$runpicasT[i]=scr1Resu[[i]][[5]]
  }
  
  results$betaEst=list()
  for(i in 1:nRef){
    results$betaEst[[i]]=scr1Resu[[i]][[6]]
  }
  
  selecList=list()
  for(i in 1:nRef){
    selecList[[i]]=scr1Resu[[i]][[1]]
  }
  
  results$yTildLongList=list()
  for(i in 1:nRef){
    results$yTildLongList[[i]]=scr1Resu[[i]][[2]]
  }
  
  selecList<- lapply(selecList, as, "sparseMatrix")
  scr1Resu=do.call(cbind, selecList)
  rm(selecList)
  
  # create count of for each predictor,each row is the count
  # of selection for a predictor
  
  countOfSelecForAPred=matrix(Matrix::rowSums(scr1Resu[1:nAlphaSelec,,drop=F]),nrow=nPredics)
  testCovCountMat=countOfSelecForAPred[testCovInd,,drop=F]
  countOfSelecForAPred=as(matrix(Matrix::colSums(testCovCountMat),nrow=1),"sparseMatrix")
  rm(scr1Resu,testCovInd)
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
