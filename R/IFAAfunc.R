##' Robust association identification and inference for absolute abundance in microbiome analyses
##'
##' Make inference on the association of covariates of microbiome
##'
##' The `IFAA()` uses a novel approach to make inference on the association of covariates with the absolute abundance (AA) of microbiome in an ecosystem.
##'
##' @param MicrobData Microbiome data matrix containing microbiome abundance with each row per sample and each column per taxon/OTU/ASV. It should contain an `"id"` variable to correspond to the `"id"` variable in the covariates data: `CovData`.
##' @param CovData Covariates data matrix containing covariates and confounders with each row per sample and each column per variable. It should also contain an `"id"` variable to correspond to the `"id"` variable in the microbiome data: `MicrobData`.
##' @param linkIDname Variable name of the `"id"` variable in both `MicrobData` and `CovData`. The two data sets will be merged by this `"id"` variable.
##' @param testCov Covariates that are of primary interest for testing and estimating the associations. It corresponds to $X_i$ in the equation. Default is `NULL` which means all covariates are `testCov`.
##' @param ctrlCov Potential confounders that will be adjusted in the model. It corresponds to $W_i$ in the equation. Default is `NULL` which means all covariates except those in `testCov` are adjusted as confounders.
##' @param testMany This takes logical value `TRUE` or `FALSE`. If `TRUE`, the `testCov` will contain all the variables in `CovData` provided `testCov` is set to be `NULL`. The default value is `TRUE` which does not do anything if `testCov` is not `NULL`.
##' @param ctrlMany This takes logical value `TRUE` or `FALSE`. If `TRUE`, all variables except `testCov` are considered as control covariates provided `ctrlCov` is set to be `NULL`. The default value is `TRUE` which does not do anything if `ctrlCov` is not `NULL`.
##' @param nRef The number of randomly picked reference taxa used in phase 1. Default number is `40`.
##' @param nPermu The number of permutation used in phase 1. Default number is `40`.
##' @param refTaxa Reference taxa will be used in phase 1. The default is `NULL` since the algorithm will pick those randomly. If a vector of taxa names are provided, the algorithm will use the provided taxa instead of randomly picking random reference taxa.
##' @param reguMethod regularization approach used in phase 1 of the algorithm. Take value `"mcp"` or `"lasso"`, default is `"mcp"`.
##' @param fwerRate The family wise error rate for identifying taxa/OTU/ASV associated with `testCov` in phase 1. Default is `0.25`.
##' @param paraJobs Number of parallel jobs that will be registered to run the algorithm. Default is `8`. If specified as `NULL`, it will automatically detect the cores to decide the number of parallel jobs.
##' @param sequentialRun This takes a logical value `TRUE` or `FALSE`. Sometimes parallel jobs can not be successfully run for unknown reasons (such as hardware issues). For example, socket related errors may pop up or some slave cores return error message instead of numerical results. In those scenarios, setting `sequentialRun = TRUE` may help, but it will take more time to run. Default is `FALSE`.
##' @param standardize This takes a logical value `TRUE` or `FALSE`. If `TRUE`, all design matrix X in phase 1 and phase 2 will be standardized in the analyses. Default is `FALSE`.
##' @param nRefMaxForEsti The maximum number of reference taxa used in phase 2. The default is `1`.
##' @param bootB Number of bootstrap samples for obtaining confidence interval of estimates in phase 2. The default is `500`.
##' @param bootLassoAlpha The significance level in phase 2. Default is `0.05`.
##' @param refReadsThresh The threshold of non-zero sequencing reads for choosing the reference taxon in phase 2. The default is `0.2` which means at least 20% non-zero sequencing reads.
##' @param SDThresh The threshold of standard deviations of sequencing reads for choosing the reference taxon in phase 2. The default is `0.5` which means the standard deviation of sequencing reads should be at least `0.5`.
##' @param balanceCut The threshold of non-zero sequencing reads in each group of a binary variable for choosing the reference taxon in phase 2. The default number is `0.2` which means at least 20% sequencing reads are non-zero in each group.
##' @param seed Random seed for reproducibility. Default is `1`.
##' @return A list containing the estimation results.
##'
##' - `analysisResults$estByCovList`: A list containing estimating results for all the variables in `testCov`. See details.
##'
##' - `covariatesData`: A dataset containing covariates and confounders used in the analyses.
##'
##' @examples
##' data(dataM)
##' dim(dataM)
##' dataM[1:5, 1:8]
##' data(dataC)
##' dim(dataC)
##' dataC[1:5, ]
##' results <- IFAA(MicrobData = dataM,
##'                 CovData = dataC,
##'                 linkIDname = "id",
##'                 testCov = c("v1", "v2"),
##'                 ctrlCov = c("v3"), nRef = 4,
##'                 nPermu = 4,
##'                 fwerRate = 0.25,
##'                 bootB = 5)
##'
##' @references \url{https://arxiv.org/pdf/1909.10101.pdf}
##' @export
##' @md





IFAA=function(
  MicrobData,
  CovData,
  linkIDname,
  testCov=NULL,
  ctrlCov=NULL,
  testMany=T,
  ctrlMany=F,
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

