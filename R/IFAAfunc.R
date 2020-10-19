##' Robust association identification and inference for absolute abundance in microbiome analyses
##'
##' Make inference on the association of covariates of microbiome
##'
##' To model the association, the following equation is used:
##' \loadmathjax
##'
##' \mjdeqn{\log(\mathcal{Y}_i^k)|\mathcal{Y}_i^k>0=\beta^{0k}+X_i^T\beta^k+W_i^T\gamma^k+Z_i^Tb_i+\epsilon_i^k,\hspace{0.2cm}k=1,...,K+1}{ASCII representation}
##' where
##' - \mjeqn{\mathcal{Y}_i^k}{} is the AA of taxa \mjeqn{k}{ASCII representation} in subject \mjeqn{i}{ASCII representation} in the entire ecosystem.
##' - \mjeqn{X_i}{} is the covariate matrix.
##' - \mjeqn{W_i}{ASCII representation} is the confounder matrix.
##' - \mjeqn{Z_i}{ASCII representation} is the design matrix for random effects.
##' - \mjeqn{\beta^k}{ASCII representation} is the regression coefficients that will be estimated and tested with the `IFAA()` function.
##'
##' The challenge in microbiome analysis is that \mjeqn{\mathcal{Y}_i^k}{ASCII representation} can not be observed. What is observed is its small proportion:
##' \mjeqn{Y_i^k=C_i\mathcal{Y}^k_i}{ascii}, where \mjeqn{C_i}{ASCII representation} is an unknown number between 0
##' and 1 that denote the observed proportion.
##'
##' The IFAA method can handle this challenge by
##' identifying and employing reference taxa. The `IFAA()` will estimate the parameter
##' \mjeqn{\beta^k}{ASCII representation} and their 95% confidence intervals. High-dimensional \mjeqn{X_i}{ASCII representation} is handled by
##' regularization.
##'
##' @param MicrobData Microbiome data matrix containing microbiome abundance with each row per
##' sample and each column per taxon/OTU/ASV. It should contain an `"id"` variable to correspond
##' to the `"id"` variable in the covariates data: `CovData`. This argument can also take file
##' directory path. For example, `MicrobData="C:\...\microbiomeData.tsv"`.
##' @param CovData Covariates data matrix containing covariates and confounders with each row
##' per sample and each column per variable. It should also contain an `"id"` variable to
##' correspond to the `"id"` variable in the microbiome data: `MicrobData`. This argument can also take file
##' directory path. For example, `CovData = "C:\...\covariatesData.tsv"`.
##' @param linkIDname Variable name of the `"id"` variable in both `MicrobData` and `CovData`. The two data sets will be merged by this `"id"` variable.
##' @param testCov Covariates that are of primary interest for testing and estimating the associations. It corresponds to $X_i$ in the equation. Default is `NULL` which means all covariates are `testCov`.
##' @param ctrlCov Potential confounders that will be adjusted in the model. It corresponds to $W_i$ in the equation. Default is `NULL` which means all covariates except those in `testCov` are adjusted as confounders.
##' @param testMany This takes logical value `TRUE` or `FALSE`. If `TRUE`, the `testCov` will contain all the variables in `CovData` provided `testCov` is set to be `NULL`. The default value is `TRUE` which does not do anything if `testCov` is not `NULL`.
##' @param ctrlMany This takes logical value `TRUE` or `FALSE`. If `TRUE`, all variables except `testCov` are considered as control covariates provided `ctrlCov` is set to be `NULL`. The default value is `TRUE` which does not do anything if `ctrlCov` is not `NULL`.
##' @param nRef The number of randomly picked reference taxa used in phase 1. Default number is `40`.
##' @param nPermu The number of permutation used in phase 1. Default number is `40`.
##' @param x1permut This takes a logical value `TRUE` or `FALSE`. If true, it will permute the variables in testCov. If false, it will use residual-permutation proposed by Freedman and Lane (1983).
##' @param refTaxa A vector of taxa or OTU or ASV names. These are reference taxa specified by the user to be used in phase 1. If the number of reference taxa is less than 'nRef', the algorithm will randomly pick extra reference taxa to make up 'nRef'. The default is `NULL` since the algorithm will pick reference taxa randomly.
##' @param reguMethod regularization approach used in phase 1 of the algorithm. Default is `"mcp"`. Other methods are under development.
##' @param fwerRate The family wise error rate for identifying taxa/OTU/ASV associated with `testCov` in phase 1. Default is `0.25`.
##' @param sequentialRun This takes a logical value `TRUE` or `FALSE`. Sometimes parallel jobs can not be successfully run for unknown reasons. For example, socket related errors may pop up or some slave cores return simple error instead of numerical results. In those scenarios, setting `sequentialRun = TRUE` may help, but it will take more time to run. Default is `FALSE`.
##' @param paraJobs If `sequentialRun` is `FALSE`, this specifies the number of parallel jobs that will be registered to run the algorithm. Default is `8`. If specified as `NULL`, it will automatically detect the cores to decide the number of parallel jobs.
##' @param standardize This takes a logical value `TRUE` or `FALSE`. If `TRUE`, all design matrix X in phase 1 and phase 2 will be standardized in the analyses. Default is `FALSE`.
##' @param nRefMaxForEsti The maximum number of reference taxa used in phase 2. The default is `1`.
##' @param allFunc all the user-defined function names that will be passed to the parallel computing environment (foreach loop).
##' @param bootB Number of bootstrap samples for obtaining confidence interval of estimates in phase 2. The default is `500`.
##' @param bootLassoAlpha The significance level in phase 2. Default is `0.05`.
##' @param refReadsThresh The threshold of non-zero sequencing reads for choosing the reference taxon in phase 2. The default is `0.2` which means at least 20% non-zero sequencing reads.
##' @param SDThresh The threshold of standard deviations of sequencing reads for choosing the reference taxon in phase 2. The default is `0.5` which means the standard deviation of sequencing reads should be at least `0.5`.
##' @param SDquantilThresh Threshold for the quantile of standard deviation for selecting final reference taxon
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
##' \dontrun{
##' results <- IFAA(MicrobData = dataM,
##'                 CovData = dataC,
##'                 linkIDname = "id",
##'                 testCov = c("v1", "v2"),
##'                 ctrlCov = c("v3"), nRef = 4,
##'                 nPermu = 4,
##'                 fwerRate = 0.25,
##'                 bootB = 5)
##'}
##'
##' @references Li et al.(2020) IFAA: Robust association identification and Inference For Absolute Abundance in microbiome analyses. arXiv:1909.10101v3
##' @references Zhang CH (2010) Nearly unbiased variable selection under minimax concave penalty. Annals of Statistics. 38(2):894-942.
##' @references Freedman and Lane (1983) A nonstochastic interpretation of reported significance levels. Journal of Business & Economic Statistics. 1(4):292-298.

##' @importFrom stats cor kmeans na.omit quantile sd
##' @importFrom utils read.csv read.table tail
##' @importFrom methods as
##' @importFrom foreach `%dopar%` foreach
##' @importFrom future availableCores
##' @importFrom Matrix Diagonal Matrix
##' @importFrom HDCI bootLOPR
##' @importFrom picasso picasso
##' @import mathjaxr
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

  if(length(refTaxa)>0){
    if(sum(refTaxa%in%microbName)!=length(refTaxa)){
      stop("Error: One or more of the specified reference taxa have no sequencing reads
      or are not in the data set. Double check the names of the reference taxa and their
           sparsity levels.")
    }
  }

  if(nRef>(length(microbName))){
    stop("Error: number of random reference taxa can not be larger than the total number
           of taxa in the data. Try lower nRef")
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

