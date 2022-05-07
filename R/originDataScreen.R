
#---------------------------------------------------------------------------------------
## Original screen function to select independent reference taxa
#---------------------------------------------------------------------------------------
originDataScreen=function(
  data,
  testCovInd,
  nRef,
  paraJobs,
  refTaxa,
  maxDimensionScr=434*5*10^5,
  sequentialRun,
  allFunc,
  Mprefix,
  covsPrefix,
  binPredInd,
  adjust_method,
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

  # overwrite nRef if the reference taxon is specified
  nRef=length(refTaxa)

  startT1=proc.time()[3]
  if(length(paraJobs)==0){
    availCores=availableCores()
    if(is.numeric(availCores))paraJobs=max(1,availableCores()-2)
    if(!is.numeric(availCores))paraJobs=1
  }
  batch=paraJobs
  forLoopN=ceiling(nRef/batch)

  if(!sequentialRun){
    message(paraJobs, " parallel jobs are registered for analyzing ", nRef, " reference taxa in Phase 1")
  }

  for (jj in 1:forLoopN){

    cl<-parallel::makeCluster(paraJobs)

    parallel::clusterExport(cl=cl, varlist=allFunc,envir=environment())

    doParallel::registerDoParallel(cl)

    if(length(seed)>0){
      set.seed(as.numeric(seed)+jj+10^4)
      parallel::clusterSetRNGStream(cl=cl,(as.numeric(seed)+jj+10^5))
    }

    if(sequentialRun){foreach::registerDoSEQ()}

    #start parallel computing
    if(forLoopN>1 & jj<forLoopN){
      scr1Resu.j=foreach(i=((jj-1)*batch+1):(jj*batch),
                         .multicombine=TRUE,
                         .packages=c("glmnet","Matrix"),
                         .errorhandling="pass") %dopar% {


                           ii=which(taxaNames==refTaxa[i])
                           dataForEst=dataRecovTrans(data=data,ref=refTaxa[i],Mprefix=Mprefix,
                                                     covsPrefix=covsPrefix)

                           xTildLongTild.i=dataForEst$xTildalong
                           yTildLongTild.i=dataForEst$UtildaLong
                           rm(dataForEst)
                           gc()

                           maxSubSamplSiz=min(50000,floor(maxDimensionScr/ncol(xTildLongTild.i)))
                           nToSamplFrom=nrow(xTildLongTild.i)

                             subSamplK=ceiling(nToSamplFrom/maxSubSamplSiz)
                             if(subSamplK==1)maxSubSamplSiz=nToSamplFrom

                             nRuns=ceiling(subSamplK/3)

                             for (k in 1:nRuns){
                               rowToKeep=sample(nToSamplFrom,maxSubSamplSiz)
                               x=as((xTildLongTild.i[rowToKeep,]),"sparseMatrix")
                               y=yTildLongTild.i[rowToKeep]


                               if (dim(x)[1]>(3*dim(x)[2])) {
                                 Penal.i=runlinear(x=x,y=y,
                                                   nPredics=nPredics)
                                 BetaNoInt.k=as((0+(Penal.i$betaNoInt!=0)),"sparseVector")
                                 EstNoInt.k<-abs(Penal.i$coef_est_noint)
                               } else {
                                 Penal.i=runGlmnet(x=x,y=y,
                                                   nPredics=nPredics)
                                 BetaNoInt.k=as((0+(Penal.i$betaNoInt!=0)),"sparseVector")
                                 EstNoInt.k<-abs(Penal.i$betaNoInt)
                               }
                               rm(Penal.i)
                               if(k==1) {
                                 BetaNoInt.i=BetaNoInt.k
                                 EstNoInt.i=EstNoInt.k
                               }
                               if(k>1) {
                                 BetaNoInt.i=BetaNoInt.i+BetaNoInt.k
                                 EstNoInt.i=EstNoInt.i+EstNoInt.k
                               }
                               rm(BetaNoInt.k,EstNoInt.k)
                             }
                             rm(k,x,y,xTildLongTild.i)
                             gc()
                             BetaNoInt.i=BetaNoInt.i/nRuns
                             EstNoInt.i=EstNoInt.i/nRuns

                           selection.i=as(rep(0,nAlphaSelec),"sparseVector")
                           coef.i<-rep(0,nAlphaSelec)
                           if (ii==1){
                             selection.i[-seq(1,nPredics)]=BetaNoInt.i
                             coef.i[-seq(1,nPredics)]=EstNoInt.i
                           }
                           if (ii==nTaxa) {
                             selection.i[-seq((nAlphaSelec-nPredics+1),nAlphaSelec)]=BetaNoInt.i
                             coef.i[-seq((nAlphaSelec-nPredics+1),nAlphaSelec)]=EstNoInt.i
                           }
                           if ((ii>1) & (ii<nTaxa)) {
                             selection.i[1:(nPredics*(ii-1))]=BetaNoInt.i[1:(nPredics*(ii-1))]
                             selection.i[(nPredics*ii+1):nAlphaSelec]=BetaNoInt.i[(nPredics*(ii-1)+1):nAlphaNoInt]
                             coef.i[1:(nPredics*(ii-1))]=EstNoInt.i[1:(nPredics*(ii-1))]
                             coef.i[(nPredics*ii+1):nAlphaSelec]=EstNoInt.i[(nPredics*(ii-1)+1):nAlphaNoInt]
                           }
                           rm(BetaNoInt.i)
                           # create return vector
                           recturnlist=list()
                           recturnlist[[1]]=selection.i
                           recturnlist[[2]]<-coef.i
                           rm(selection.i,yTildLongTild.i,coef.i)
                           return(recturnlist)
                         }
      parallel::stopCluster(cl)
      gc()
      if(jj==1)scr1Resu=scr1Resu.j
      if(jj>1)scr1Resu=do.call(base::c,list(scr1Resu,scr1Resu.j))
    }

    if(jj==forLoopN){
      scr1Resu.j=foreach(i=((forLoopN-1)*batch+1):nRef,
                         .multicombine=TRUE,
                         .packages=c("glmnet","Matrix"),
                         .errorhandling="pass") %dopar% {


                           ii=which(taxaNames==refTaxa[i])
                           dataForEst=dataRecovTrans(data=data,ref=refTaxa[i],Mprefix=Mprefix,
                                                     covsPrefix=covsPrefix)

                           xTildLongTild.i=dataForEst$xTildalong
                           yTildLongTild.i=dataForEst$UtildaLong
                           rm(dataForEst)
                           gc()

                           maxSubSamplSiz=min(50000,floor(maxDimensionScr/ncol(xTildLongTild.i)))
                           nToSamplFrom=nrow(xTildLongTild.i)

                             subSamplK=ceiling(nToSamplFrom/maxSubSamplSiz)
                             if(subSamplK==1)maxSubSamplSiz=nToSamplFrom

                             nRuns=ceiling(subSamplK/3)

                             for (k in 1:nRuns){
                               rowToKeep=sample(nToSamplFrom,maxSubSamplSiz)
                               x=as((xTildLongTild.i[rowToKeep,]),"sparseMatrix")
                               y=(yTildLongTild.i[rowToKeep])

                               if (dim(x)[1]>(3*dim(x)[2])) {
                                 Penal.i=runlinear(x=x,y=y,
                                                   nPredics=nPredics)
                                 BetaNoInt.k=as((0+(Penal.i$betaNoInt!=0)),"sparseVector")
                                 EstNoInt.k<-abs(Penal.i$coef_est_noint)
                               } else {
                                 Penal.i=runGlmnet(x=x,y=y,
                                                   nPredics=nPredics)
                                 BetaNoInt.k=as((0+(Penal.i$betaNoInt!=0)),"sparseVector")
                                 EstNoInt.k<-abs(Penal.i$betaNoInt)
                               }

                               rm(Penal.i)
                               gc()
                               if(k==1) {
                                 BetaNoInt.i=BetaNoInt.k
                                 EstNoInt.i=EstNoInt.k
                               }
                               if(k>1) {
                                 BetaNoInt.i=BetaNoInt.i+BetaNoInt.k
                                 EstNoInt.i=EstNoInt.i+EstNoInt.k
                               }
                               rm(BetaNoInt.k,EstNoInt.k)
                             }
                             rm(k,x,y,xTildLongTild.i)
                             gc()
                             BetaNoInt.i=BetaNoInt.i/nRuns
                             EstNoInt.i=EstNoInt.i/nRuns

                           selection.i=as(rep(0,nAlphaSelec),"sparseVector")
                           coef.i<-rep(0,nAlphaSelec)
                           if (ii==1){
                             selection.i[-seq(1,nPredics)]=BetaNoInt.i
                             coef.i[-seq(1,nPredics)]=EstNoInt.i
                           }
                           if (ii==nTaxa) {
                             selection.i[-seq((nAlphaSelec-nPredics+1),nAlphaSelec)]=BetaNoInt.i
                             coef.i[-seq((nAlphaSelec-nPredics+1),nAlphaSelec)]=EstNoInt.i
                           }
                           if ((ii>1) & (ii<nTaxa)) {
                             selection.i[1:(nPredics*(ii-1))]=BetaNoInt.i[1:(nPredics*(ii-1))]
                             selection.i[(nPredics*ii+1):nAlphaSelec]=BetaNoInt.i[(nPredics*(ii-1)+1):nAlphaNoInt]
                             coef.i[1:(nPredics*(ii-1))]=EstNoInt.i[1:(nPredics*(ii-1))]
                             coef.i[(nPredics*ii+1):nAlphaSelec]=EstNoInt.i[(nPredics*(ii-1)+1):nAlphaNoInt]
                           }
                           rm(BetaNoInt.i)
                           # create return vector
                           recturnlist=list()
                           recturnlist[[1]]=selection.i
                           recturnlist[[2]]<-coef.i
                           rm(selection.i,yTildLongTild.i,coef.i)
                           return(recturnlist)
                         }
      parallel::stopCluster(cl)
      if(forLoopN==1)scr1Resu=scr1Resu.j
      if(forLoopN>1)scr1Resu=do.call(base::c,list(scr1Resu,scr1Resu.j))
      gc()
    }

    # if(jj>0 & (jj%%(ceiling(forLoopN/10))==0)){
    #   message(round(100*jj/forLoopN,0), " percent of phase 1 analysis has been done")
    # }
  }
  rm(data)
  endT=proc.time()[3]


  selecList=list()
  for(i in 1:nRef){
    selecList[[i]]=scr1Resu[[i]][[1]]
  }




  estList=list()
  for(i in 1:nRef){
    estList[[i]]=scr1Resu[[i]][[2]]
  }

  rm(scr1Resu)

  selecList<- lapply(selecList, as, "sparseMatrix")
  estList<-lapply(estList,as.matrix)
  scr1ResuSelec=do.call(cbind, selecList)
  scr1ResuEst<-do.call(cbind,estList)

  rm(selecList,estList)



  # create count of selection for individual testCov
  countOfSelecForAllPred=as(matrix(Matrix::rowSums(scr1ResuSelec),nrow=nPredics),"sparseMatrix")
  EstOfAllPred<-matrix(rowSums(scr1ResuEst),nrow=nPredics)

  testCovCountMat=countOfSelecForAllPred[testCovInd,,drop=FALSE]
  testEstMat<-EstOfAllPred[testCovInd,,drop=FALSE]
  rm(scr1ResuSelec,testCovInd,countOfSelecForAllPred,EstOfAllPred)

  # create overall count of selection for all testCov as a whole
  countOfSelecForAPred=as(matrix(Matrix::colSums(testCovCountMat),nrow=1),"sparseMatrix")
  estOfSelectForAPred<-matrix(Matrix::colSums(testEstMat),nrow=1)
  gc()

  colnames(countOfSelecForAPred)=taxaNames
  colnames(estOfSelectForAPred)<-taxaNames
  rm(taxaNames)

  # return results
  results$testCovCountMat=testCovCountMat
  results$testEstMat<-testEstMat
  rm(testCovCountMat,testEstMat)
  results$countOfSelecForAPred=countOfSelecForAPred
  results$estOfSelectForAPred<-estOfSelectForAPred
  rm(countOfSelecForAPred,estOfSelectForAPred)
  return(results)
}
