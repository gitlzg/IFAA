##' @export


XbetaAndResidu=function(
  method,
  data,
  testCovInd,
  nRef,
  paraJobs,
  lambda=NULL,
  refTaxa,
  sequentialRun,
  standardize,
  allFunc,
  Mprefix,
  covsPrefix,
  binPredInd,
  seed){
  
  results=list()
  
  # generate x matrix for the permutation
  dataForEst1=dataRecovTrans(data=data,ref=refTaxa[1],Mprefix=Mprefix,
                             covsPrefix=covsPrefix)
  results$xTildLong=dataForEst1$xTildalong
  rm(dataForEst1)
  
  # generate reduced data
  basicInfo=dataInfo(data=data,Mprefix=Mprefix,
                     covsPrefix=covsPrefix,
                     binPredInd=binPredInd)
  
  predNames=basicInfo$predNames
  testCovNam=predNames[testCovInd]
  colnames(data)%in%testCovNam
  reducedData=data[,!(colnames(data)%in%testCovNam)]
  rm(basicInfo)
  
  # load reduced data info
  reducedBasicInfo=dataInfo(data=reducedData,Mprefix=Mprefix,
                            covsPrefix=covsPrefix,
                            binPredInd=binPredInd)
  
  taxaNames=reducedBasicInfo$taxaNames
  nTaxa=reducedBasicInfo$nTaxa
  nPredics=reducedBasicInfo$nPredics
  rm(reducedBasicInfo)
  gc()
  
  nNorm=nTaxa-1
  nAlphaNoInt=nPredics*nNorm
  nAlphaSelec=nPredics*nTaxa
  
  countOfSelec=rep(0,nAlphaSelec)
  resultsByRefTaxon=list()
  
  # overwrite nRef if the reference taxon is specified
  nRef=length(refTaxa)
  
  startT=proc.time()[3]
  cat("start generating residues for permutation","\n")
  if(length(paraJobs)==0){
    availCores=availableCores()
    if(is.numeric(availCores))paraJobs=max(1,availableCores()-2)
    if(!is.numeric(availCores))paraJobs=1
  }
  
  c4<-snow::makeCluster(paraJobs)
  
  if(!sequentialRun){
    cat(paraJobs, "parallel jobs are registered for generate residues in Phase 1a.","\n")
  }
  
  snow::clusterExport(c4, allFunc)
  doSNOW::registerDoSNOW(c4)
  
  if(sequentialRun){foreach::registerDoSEQ()}
  
  startT1=proc.time()[3]
  # start parallel computing
  residu1Resu=foreach(i=1:nRef,.multicombine=T,
                      .packages=c("picasso","expm","doSNOW","snow","foreach","Matrix"),
                      .errorhandling="pass") %dopar% {
                        
                        ii=which(taxaNames==refTaxa[i])
                        dataForEst=dataRecovTrans(data=reducedData,ref=refTaxa[i],Mprefix=Mprefix,
                                                  covsPrefix=covsPrefix)
                        xTildLongTild.i=dataForEst$xTildalong
                        yTildLongTild.i=dataForEst$UtildaLong
                        rm(dataForEst)
                        gc()
                        
                        if(method=="mcp") {
                          Penal.i=runPicasso(x=xTildLongTild.i,y=yTildLongTild.i,
                                             lambda=lambda,nPredics=nPredics,
                                             standardize=standardize,
                                             method="mcp",permutY=F,seed=seed,seedi=i)
                        }
                        
                        betaInt=Penal.i$betaInt
                        overalIntercp=Penal.i$overalIntercp
                        rm(Penal.i)
                        gc()
                        
                        xBeta=xTildLongTild.i%*%betaInt+overalIntercp
                        rm(xTildLongTild.i)
                        
                        residu=yTildLongTild.i-xBeta
                        rm(yTildLongTild.i)
                        
                        # create return vector
                        recturnlist=list()
                        recturnlist[[1]]=xBeta
                        recturnlist[[2]]=residu
                        return(recturnlist)
                      }
  snow::stopCluster(c4)
  rm(data)
  gc()
  
  endT=proc.time()[3]
  
  cat("Generating residu is done and took",(endT-startT1)/60,"minutes","\n")
  
  xBetaList=list()
  for(i in 1:nRef){
    xBetaList[[i]]=residu1Resu[[i]][[1]]
  }
  
  residuList=list()
  for(i in 1:nRef){
    residuList[[i]]=residu1Resu[[i]][[2]]
  }
  rm(residu1Resu)
  
  # return results
  results$xBetaList=xBetaList
  rm(xBetaList)
  results$residuList=residuList
  rm(residuList)
  return(results)
}
