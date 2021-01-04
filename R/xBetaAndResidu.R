#---------------------------------------------------------------------------------------
## Function to calculate "X\beta" and residuals for reduced models
## to prepare for Freedman and Lane (1983) permutation
#---------------------------------------------------------------------------------------

XbetaAndResidu=function(
  method,
  data,
  testCovInd,
  nRef,
  paraJobs,
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
  message("start generating residues for permutation")
  if(length(paraJobs)==0){
    availCores=availableCores()
    if(is.numeric(availCores))paraJobs=max(1,availableCores()-2)
    if(!is.numeric(availCores))paraJobs=1
  }
  if(!sequentialRun){
    message(paraJobs, "parallel jobs are registered for generate residues in Phase 1a.")
  }
  batch=paraJobs
  forLoopN=ceiling(nRef/batch)

  for(jj in 1:forLoopN){
    cl<-parallel::makeCluster(paraJobs)
    
    parallel::clusterExport(cl=cl, varlist=allFunc,envir=parent.env(environment()))
    doParallel::registerDoParallel(cl)
    
    if(length(seed)>0){
      set.seed(as.numeric(seed)+jj+10^8)
      parallel::clusterSetRNGStream(cl=cl,(as.numeric(seed)+jj+10^9))
    }
    if(sequentialRun){foreach::registerDoSEQ()}
    
    startT1=proc.time()[3]
    if(forLoopN>1 & jj<forLoopN){
      # start parallel computing
      residu1Resu.j=foreach(i=((jj-1)*batch+1):(jj*batch),.multicombine=T,
                            .packages=c("picasso","Matrix"),
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
                                                   nPredics=nPredics,
                                                   method="mcp",permutY=FALSE,
                                                   standardize=standardize,
                                                   seed=seed,seedi=i)
                                
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
      parallel::stopCluster(cl)
      rm(data)
      gc()
      if(jj==1)residu1Resu=residu1Resu.j
      if(jj>1)residu1Resu=do.call(c,list(residu1Resu,residu1Resu.j))
    }
    if(jj==forLoopN){
      # start parallel computing
      residu1Resu.j=foreach(i=((jj-1)*batch+1):nRef,.multicombine=T,
                            .packages=c("picasso","Matrix"),
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
                                                   nPredics=nPredics,
                                                   method="mcp",permutY=FALSE,
                                                   standardize=standardize,
                                                   seed=seed,seedi=i)

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
      parallel::stopCluster(cl)
      rm(data)
      gc()
      if(forLoopN==1)residu1Resu=residu1Resu.j
      if(forLoopN>1)residu1Resu=do.call(c,list(residu1Resu,residu1Resu.j))
    }
  }
  endT=proc.time()[3]
  
  message("Generating residu is done and took ",round((endT-startT1)/60,3),"minutes")
  
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
