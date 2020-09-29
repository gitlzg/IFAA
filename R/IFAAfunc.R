
##' @export




IFAA=function(
  MicrobData,
  CovData,
  linkIDname,
  testCov=NULL,
  ctrlCov=NULL,
  testMany=T,
  ctrlMany=T,
  nRef=40,
  nRefMaxForEsti=1,
  nPermu=40,
  x1permut=T,
  refTaxa=NULL,
  reguMethod=c("mcp"), 
  fwerRate=0.25,
  paraJobs=NULL,
  bootB=500,
  bootLassoAlpha=0.05,
  standardize=F,
  sequentialRun=F,
  allFunc=allUserFunc(),
  refReadsThresh=0.2,
  SDThresh=0.05,
  SDquantilThresh=0,
  balanceCut=0.2,
  seed=1
){
  results=list()
  start.time = proc.time()[3] 
  runMeta=metaData(MicrobData=MicrobData,CovData=CovData,
                   linkIDname=linkIDname,testCov=testCov,
                   ctrlCov=ctrlCov,testMany=testMany,ctrlMany=ctrlMany)
  data=runMeta$data
  results$covariatesData=runMeta$covariatesData
  binaryInd=runMeta$binaryInd
  covsPrefix=runMeta$covsPrefix
  Mprefix=runMeta$Mprefix
  testCovInd=runMeta$testCovInd
  testCovInOrder=runMeta$testCovInOrder
  testCovInNewNam=runMeta$testCovInNewNam
  ctrlCov=runMeta$ctrlCov
  microbName=runMeta$microbName
  results$covriateNames=runMeta$xNames
  rm(runMeta)
  
  # if(length(refTaxa)>1){
  #   stop("Error: Too many reference taxa specified. Only one at a time is allowed.")
  #   }
  if(length(refTaxa)>0){
    nRef=length(refTaxa)
    if(sum(refTaxa%in%microbName)!=length(refTaxa)){
      stop("Error: One or more of the specified reference taxa is not in the data set, 
         double check the name of the reference taxon.")
    }
  }
  results$analysisResults=Regulariz(data=data,testCovInd=testCovInd,
                                    testCovInOrder=testCovInOrder,
                                    testCovInNewNam=testCovInNewNam,
                                    microbName=microbName,nRef=nRef,
                                    nRefMaxForEsti=nRefMaxForEsti,
                                    nPermu=nPermu,binaryInd=binaryInd,
                                    x1permut=x1permut,
                                    covsPrefix=covsPrefix,Mprefix=Mprefix,
                                    refTaxa=refTaxa,paraJobs=paraJobs,
                                    reguMethod=reguMethod,fwerRate=fwerRate,
                                    bootB=bootB,bootLassoAlpha=bootLassoAlpha,
                                    standardize=standardize,
                                    sequentialRun=sequentialRun,
                                    allFunc=allFunc,refReadsThresh=refReadsThresh,
                                    SDThresh=SDThresh,
                                    SDquantilThresh=SDquantilThresh,
                                    balanceCut=balanceCut,seed=seed
  )
  rm(data)
  
  results$testCov=testCovInOrder
  results$ctrlCov=ctrlCov
  results$microbName=microbName
  results$bootB=bootB
  results$bootLassoAlpha=bootLassoAlpha
  results$refReadsThresh=refReadsThresh
  results$balanceCut=balanceCut
  results$SDThresh=SDThresh
  results$paraJobs=paraJobs
  results$SDquantilThresh=SDquantilThresh
  results$nRef=nRef
  results$nPermu=nPermu
  results$x1permut=x1permut
  
  if(length(seed)==1){
    results$seed=seed
  }else{
    results$seed="No seed used."
  }
  
  rm(testCovInOrder,ctrlCov,microbName)
  
  totalTimeMins = (proc.time()[3] - start.time)/60
  cat("The entire analysis took",totalTimeMins, "minutes","\n")
  
  results$totalTimeMins=totalTimeMins
  
  return(results)
}

