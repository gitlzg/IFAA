##' Robust association identification and inference for absolute abundance in microbiome analyses
##'
##' Make inference on the association of microbiome with covariates
##'
##' Most of the time, users just need to feed the first five inputs to the function: `MicrobData`, `CovData`, `linkIDname`, `testCov` and `ctrlCov`. All other inputs can just take their default values. 
##' To model the association, the following equation is used:
##'
##' \loadmathjax
##'
#' \mjdeqn{\log(\mathcal{Y}_i^k)|\mathcal{Y}_i^k>0=\beta^{0k}+X_i^T\beta^k+W_i^T\gamma^k+Z_i^Tb_i+\epsilon_i^k,\hspace{0.2cm}k=1,...,K+1}{}
##' where
##' - \mjeqn{\mathcal{Y}_i^k}{} is the AA of taxa \mjeqn{k}{} in subject \mjeqn{i}{} in the entire ecosystem.
##' - \mjeqn{X_i}{} is the covariate matrix.
##' - \mjeqn{W_i}{} is the confounder matrix.
##' - \mjeqn{Z_i}{} is the design matrix for random effects.
##' - \mjeqn{\beta^k}{} is the regression coefficients that will be estimated and tested with the `IFAA()` function.
##'
##' The challenge in microbiome analysis is that \mjeqn{\mathcal{Y}_i^k}{} can not be observed. What is observed is its small proportion:
##' \mjeqn{Y_i^k=C_i\mathcal{Y}^k_i}{}, where \mjeqn{C_i}{} is an unknown number between 0
##' and 1 that denote the observed proportion.
##'
##' The IFAA method can successfully addressed this challenge. The `IFAA()` will estimate the parameter
##' \mjeqn{\beta^k}{} and their 95% confidence intervals. High-dimensional \mjeqn{X_i}{} is handled by
##' regularization.
##'
##' @param MicrobData Microbiome data matrix containing microbiome abundance with each row per
##' sample and each column per taxon/OTU/ASV. It should contain an `"id"` variable to correspond
##' to the `"id"` variable in the covariates data: `CovData`. This argument can take
##' directory path. For example, `MicrobData="C://...//microbiomeData.tsv"`.
##' @param CovData Covariates data matrix containing covariates and confounders with each row
##' per sample and each column per variable. It should also contain an `"id"` variable to
##' correspond to the `"id"` variable in the microbiome data: `MicrobData`. This argument can take
##' directory path. For example, `CovData = "C://...//covariatesData.tsv"`.
##' @param linkIDname Variable name of the `"id"` variable in both `MicrobData` and `CovData`. The two data sets will be merged by this `"id"` variable.
##' @param testCov Covariates that are of primary interest for testing and estimating the associations. It corresponds to $X_i$ in the equation. Default is `NULL` which means all covariates are `testCov`.
##' @param ctrlCov Potential confounders that will be adjusted in the model. It corresponds to $W_i$ in the equation. Default is `NULL` which means all covariates except those in `testCov` are adjusted as confounders.
##' @param testMany This takes logical value `TRUE` or `FALSE`. If `TRUE`, the `testCov` will contain all the variables in `CovData` provided `testCov` is set to be `NULL`. The default value is `TRUE` which does not do anything if `testCov` is not `NULL`.
##' @param ctrlMany This takes logical value `TRUE` or `FALSE`. If `TRUE`, all variables except `testCov` are considered as control covariates provided `ctrlCov` is set to be `NULL`. The default value is `FALSE`.
##' @param nRef The number of randomly picked reference taxa used in phase 1. Default number is `40`.
##' @param nRefMaxForEsti The maximum number of final reference taxa used in phase 2. The default is `2`.
##' @param refTaxa A vector of taxa or OTU or ASV names. These are reference taxa specified by the user to be used in phase 1. If the number of reference taxa is less than 'nRef', the algorithm will randomly pick extra reference taxa to make up 'nRef'. The default is `NULL` since the algorithm will pick reference taxa randomly.
##' @param adjust_method The adjusting method for p value adjustment. Default is "BY" for dependent FDR adjustment. It can take any adjustment method for p.adjust function in R.
##' @param fdrRate The false discovery rate for identifying taxa/OTU/ASV associated with `testCov`. Default is `0.15`.
##' @param paraJobs If `sequentialRun` is `FALSE`, this specifies the number of parallel jobs that will be registered to run the algorithm. If specified as `NULL`, it will automatically detect the cores to decide the number of parallel jobs. Default is `NULL`. 
##' @param bootB Number of bootstrap samples for obtaining confidence interval of estimates in phase 2 for the high dimensional regression. The default is `500`.
##' @param standardize This takes a logical value `TRUE` or `FALSE`. If `TRUE`, all design matrix X in phase 1 and phase 2 will be standardized in the analyses. Default is `FALSE`.
##' @param sequentialRun This takes a logical value `TRUE` or `FALSE`. Default is `FALSE`. This argument could be useful for debug.
##' @param refReadsThresh The threshold of proportion of non-zero sequencing reads for choosing the reference taxon in phase 2. The default is `0.2` which means at least 20% non-zero sequencing reads.
##' @param taxkeepThresh The threshold of number of non-zero sequencing reads for each taxon to be included into the analysis. The default is `0` which means taxon with at least `0` sequencing reads will be included into the analysis
##' @param SDThresh The threshold of standard deviations of sequencing reads for been chosen as the reference taxon in phase 2. The default is `0.05` which means the standard deviation of sequencing reads should be at least `0.05` in order to be chosen as reference taxon.
##' @param SDquantilThresh The threshold of the quantile of standard deviation of sequencing reads, above which could be selected as reference taxon. The default is `0`.
##' @param balanceCut The threshold of the proportion of non-zero sequencing reads in each group of a binary variable for choosing the final reference taxa in phase 2. The default number is `0.2` which means at least 20% non-zero sequencing reads in each group are needed to be eligible for being chosen as a final reference taxon.
##' @param seed Random seed for reproducibility. Default is `1`. It can be set to be NULL to remove seeding.
##' @return A list containing the estimation results.
##'
##' - `sig_results`: A list containing estimating results that are statistically significant.
##' - `full_results`: A list containing all estimating results. NA denotes unestimable.
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
##' \donttest{
##' results <- IFAA(MicrobData = dataM,
##'                 CovData = dataC,
##'                 linkIDname = "id",
##'                 testCov = c("v1", "v2"),
##'                 ctrlCov = c("v3"), 
##'                 fdrRate = 0.15)
##'
##'}
##'
##'
##'
##' @references Li et al.(2021) IFAA: Robust association identification and Inference For Absolute Abundance in microbiome analyses. Journal of the American Statistical Association
##' @references Zhang CH (2010) Nearly unbiased variable selection under minimax concave penalty. Annals of Statistics. 38(2):894-942.
##' @references Liu et al.(2020) A bootstrap lasso + partial ridge method to construct confidence intervals for parameters in high-dimensional sparse linear models. Statistica Sinica
##' @importFrom methods as
##' @importFrom foreach `%dopar%` foreach
##' @importFrom future availableCores
##' @importFrom Matrix Diagonal Matrix
##' @importFrom HDCI bootLOPR
##' @importFrom qlcMatrix corSparse
##' @import expm
##' @import rlecuyer
##' @import mathjaxr
##' @import glmnet
##' @import stats
##' @import utils
##' @export
##' @md

## arg deleted: nperm, x1permt, reguMethod, method, bootLassoAlpha


IFAA=function(
  MicrobData,
  CovData,
  linkIDname,
  testCov=NULL,
  ctrlCov=NULL,
  testMany=TRUE,
  ctrlMany=FALSE,
  nRef=40,
  nRefMaxForEsti=2,
  refTaxa=NULL,
  adjust_method="BY",
  fdrRate=0.15,
  paraJobs=NULL,
  bootB=500,
  standardize=FALSE,
  sequentialRun=FALSE,
  refReadsThresh=0.2,
  taxkeepThresh=1,
  SDThresh=0.05,
  SDquantilThresh=0,
  balanceCut=0.2,
  seed=1
){
  allFunc=allUserFunc()

  results=list()
  start.time = proc.time()[3]
  runMeta=metaData(MicrobData=MicrobData,CovData=CovData,
                   linkIDname=linkIDname,testCov=testCov,
                   ctrlCov=ctrlCov,testMany=testMany,
                   ctrlMany=ctrlMany,taxkeepThresh=taxkeepThresh)
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
  newMicrobNames=runMeta$newMicrobNames
  results$covriateNames=runMeta$xNames
  rm(runMeta)

  if(length(refTaxa)>0){
    if(sum(refTaxa%in%microbName)!=length(refTaxa)){
      stop("Error: One or more of the specified reference taxa in phase 1 have no sequencing reads
      or are not in the data set. Double check the names of the reference taxa and their
           sparsity levels.")
    }
  }


  if (nRefMaxForEsti<2) {
    nRefMaxForEsti<-2
    warning("Warning: Needs at least two final reference taxon for estimation.")
  }


  if(nRef>(length(microbName))){
    stop("Error: number of random reference taxa can not be larger than the total number
           of taxa in the data. Try lower nRef")
  }

  refTaxa_newNam=newMicrobNames[microbName%in%refTaxa]


  results$analysisResults=Regulariz(data=data,testCovInd=testCovInd,
                                    testCovInOrder=testCovInOrder,
                                    testCovInNewNam=testCovInNewNam,
                                    microbName=microbName,nRef=nRef,
                                    nRefMaxForEsti=nRefMaxForEsti,
                                    binaryInd=binaryInd,
                                    covsPrefix=covsPrefix,Mprefix=Mprefix,
                                    refTaxa=refTaxa_newNam,
                                    paraJobs=paraJobs,
                                    adjust_method=adjust_method,
                                    fwerRate=fdrRate,
                                    bootB=bootB,
                                    standardize=standardize,
                                    sequentialRun=sequentialRun,
                                    allFunc=allFunc,refReadsThresh=refReadsThresh,
                                    SDThresh=SDThresh,
                                    SDquantilThresh=SDquantilThresh,
                                    balanceCut=balanceCut,seed=seed
  )
  rm(data)

  results$sig_results<-results$analysisResults$sig_results
  results$full_results<-results$analysisResults$full_results

  results$testCov=testCovInOrder
  results$ctrlCov=ctrlCov
  results$microbName=microbName
  results$bootB=bootB
  results$refReadsThresh=refReadsThresh
  results$balanceCut=balanceCut
  results$SDThresh=SDThresh
  results$paraJobs=paraJobs
  results$SDquantilThresh=SDquantilThresh
  results$nRef=nRef

  if(length(seed)==1){
    results$seed=seed
  }else{
    results$seed="No seed used."
  }

  rm(testCovInOrder,ctrlCov,microbName)

  totalTimeMins = (proc.time()[3] - start.time)/60
  message("The entire analysis took ",round(totalTimeMins,2), " minutes")

  results$totalTimeMins=totalTimeMins

  return(results)
}


