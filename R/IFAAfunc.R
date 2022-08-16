##' Robust association identification and inference for absolute abundance in microbiome analyses
##'
##' @description
##' The IFAA function is to make inference on the association of microbiome with covariates
##'
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
##' When using this function, most of the time, users just need to feed the first three inputs to the function: `experiment_dat`, `testCov` and `ctrlCov`. All other inputs can just take their default values.
##'
##' @param experiment_dat A SummarizedExperiment object containing microbiome data and covariates (see example on how to create a SummarizedExperiment object). The microbiome data can be absolute abundance or relative abundance
##' with each column per sample and each row per taxon/OTU/ASV (or any other unit). No imputation is needed for zero-valued data points. The covariates data contains covariates and confounders with each row per sample and each
##' column per variable. The covariates data has to be numeric or binary. Categorical variables should be converted into dummy variables.
##' @param testCov Covariates that are of primary interest for testing and estimating the associations. It corresponds to $X_i$ in the equation. Default is `NULL` which means all covariates are `testCov`.
##' @param ctrlCov Potential confounders that will be adjusted in the model. It corresponds to $W_i$ in the equation. Default is `NULL` which means all covariates except those in `testCov` are adjusted as confounders.
##' @param sampleIDname Name of the sample ID variable in the data. In the case that the data does not have an ID variable, this can be ignored. Default is NULL.
##' @param testMany This takes logical value `TRUE` or `FALSE`. If `TRUE`, the `testCov` will contain all the variables in `CovData` provided `testCov` is set to be `NULL`. The default value is `TRUE` which does not do anything if `testCov` is not `NULL`.
##' @param ctrlMany This takes logical value `TRUE` or `FALSE`. If `TRUE`, all variables except `testCov` are considered as control covariates provided `ctrlCov` is set to be `NULL`. The default value is `FALSE`.
##' @param nRef The number of randomly picked reference taxa used in phase 1. Default number is `40`.
##' @param nRefMaxForEsti The maximum number of final reference taxa used in phase 2. The default is `2`.
##' @param refTaxa A vector of taxa or OTU or ASV names. These are reference taxa specified by the user to be used in phase 1. If the number of reference taxa is less than 'nRef', the algorithm will randomly pick extra reference taxa to make up 'nRef'. The default is `NULL` since the algorithm will pick reference taxa randomly.
##' @param adjust_method The adjusting method for p value adjustment. Default is "BY" for dependent FDR adjustment. It can take any adjustment method for p.adjust function in R.
##' @param fdrRate The false discovery rate for identifying taxa/OTU/ASV associated with `testCov`. Default is `0.05`.
##' @param paraJobs If `sequentialRun` is `FALSE`, this specifies the number of parallel jobs that will be registered to run the algorithm. If specified as `NULL`, it will automatically detect the cores to decide the number of parallel jobs. Default is `NULL`.
##' @param bootB Number of bootstrap samples for obtaining confidence interval of estimates in phase 2 for the high dimensional regression. The default is `500`.
##' @param standardize This takes a logical value `TRUE` or `FALSE`. If `TRUE`, the design matrix for X will be standardized in the analyses and the results. Default is `FALSE`.
##' @param sequentialRun This takes a logical value `TRUE` or `FALSE`. Default is `FALSE`. This argument could be useful for debug.
##' @param refReadsThresh The threshold of proportion of non-zero sequencing reads for choosing the reference taxon in phase 2. The default is `0.2` which means at least 20% non-zero sequencing reads.
##' @param taxDropThresh The threshold of number of non-zero sequencing reads for each taxon to be dropped from the analysis. The default is `0` which means taxon without any sequencing reads will be dropped from the analysis.
##' @param SDThresh The threshold of standard deviations of sequencing reads for been chosen as the reference taxon in phase 2. The default is `0.05` which means the standard deviation of sequencing reads should be at least `0.05` in order to be chosen as reference taxon.
##' @param SDquantilThresh The threshold of the quantile of standard deviation of sequencing reads, above which could be selected as reference taxon. The default is `0`.
##' @param balanceCut The threshold of the proportion of non-zero sequencing reads in each group of a binary variable for choosing the final reference taxa in phase 2. The default number is `0.2` which means at least 20% non-zero sequencing reads in each group are needed to be eligible for being chosen as a final reference taxon.
##' @param verbose Whether the process message is printed out to the console. The default is TRUE.
##' @return A list containing 2 elements
##' \itemize{
##' \item {`full_results`: The main results for IFAA containing the estimation and testing results for all associations between all taxa and all test covariates in `testCov`. It is a dataframe with each row
##' representing an association, and eight columns named as "taxon", "cov", "estimate", "SE.est", "CI.low", "CI.up", "adj.p.value", and "sig_ind". The columns correspond to taxon name, covariate name, association estimates,
##' standard error estimates, lower bound and upper bound of the 95% confidence interval, adjusted p value, and the indicator showing whether the association is significant after multiple testing adjustment.
##' }
##' \item {`metadata`: The metadata is a list.
##' \itemize{
##' \item {`covariatesData`: A dataset containing covariates and confounders used in the analyses.}
##' \item {`final_ref_taxon`: The final 2 reference taxon used for analysis.}
##' \item {`ref_taxon_count`: The counts of selection for the associations of all taxa with test covariates in Phase 1.}
##' \item {`totalTimeMins`: The average magnitude estimates for the associations of all taxa with test covariates in Phase 1.}
##' \item {`ref_taxon_est`: Total time used for the entire analysis.}
##' \item {`fdrRate`: FDR rate used for the analysis.}
##' \item {`adjust_method`: Multiple testing adjust method used for the analysis.}

##' }
##' }
##' }
##' @examples
##'
##' library(IFAA)
##' ## A makeup example data from Scratch. 60 taxon, 40 subjects, 3 covariates
##'
##' set.seed(1)
##' 
##' ## create an ID variable for the example data
##' ID=seq_len(40)
##' 
##' ## generate three covariates x1, x2, and x3, with x2 binary
##' x1<-rnorm(40)
##' x2<-rbinom(40,1,0.5)
##' x3<-rnorm(40)
##' dataC<-data.frame(cbind(ID,x1,x2,x3))
##'
##' ## Coefficients for x1, x2, and x3 among 60 taxa.
##' beta_1<-c(0.1,rep(0,59))
##' beta_2<-c(0,0.2,rep(0,58))
##' beta_3<-rnorm(60)
##' beta_mat<-cbind(beta_1,beta_2,beta_3)
##'
##' ## Generate absolute abundance for 60 taxa in ecosystem.
##' dataM_eco<-floor(exp(10+as.matrix(dataC[,-1])%*%t(beta_mat) + rnorm(2400,sd=0.05)))
##'
##' ## Generate sequence depth and generate observed abundance
##' Ci<-runif(40,0.01,0.05)
##' dataM<-floor(apply(dataM_eco,2,function(x) x*Ci))
##' colnames(dataM)<-paste0("rawCount",1:60)
##'
##' ## Randomly introduce 0 to make 25% sparsity level.
##' dataM[sample(seq_len(length(dataM)),length(dataM)/4)]<-0
##'
##' dataM<-data.frame(cbind(ID,dataM))
##' 
##' ## The following steps are to create a SummarizedExperiment object.
##' ## If you already have a SummarizedExperiment format data, you can
##' ## ignore the following steps and directly feed it to the IFAA function.
##' 
##' ## Merge two dataset by ID variable
##' data_merged<-merge(dataM,dataC,by="ID",all=FALSE)
##' 
##' ## Seperate microbiome data and covariate data, drop ID variable from microbiome data
##' dataM_sub<-data_merged[,colnames(dataM)[!colnames(dataM)%in%c("ID")]]
##' dataC_sub<-data_merged[,colnames(dataC)]
##' 
##' ## Create SummarizedExperiment object
##' test_dat<-SummarizedExperiment::SummarizedExperiment(
##' assays=list(MicrobData=t(dataM_sub)), colData=dataC_sub)
##'
##' ## Again, if you already have a SummarizedExperiment format data, you can
##' ## ignore the above steps and directly feed it to the IFAA function.
##'
##' set.seed(123) # For full reproducibility
##'
##' results <- IFAA(experiment_dat = test_dat,
##'                 testCov = c("x1", "x2"),
##'                 ctrlCov = c("x3"),
##'                 sampleIDname="ID",
##'                 fdrRate = 0.05)
##' ## to extract all results:
##' summary_res<-results$full_results
##' ## to extract significant results:
##' sig_results=subset(summary_res,sig_ind==TRUE)
##'
##'
##'
##' @references Li et al.(2021) IFAA: Robust association identification and Inference For Absolute Abundance in microbiome analyses. Journal of the American Statistical Association. 116(536):1595-1608

##' @importFrom foreach foreach %dopar% registerDoSEQ
##' @importFrom doRNG %dorng%
##' @importFrom parallelly availableCores
##' @importFrom parallel makeCluster clusterExport stopCluster clusterSetRNGStream
##' @importFrom doParallel registerDoParallel
##' @importFrom Matrix Diagonal Matrix sparseVector
##' @importFrom glmnet glmnet
##' @importFrom HDCI bootLOPR
##' @import mathjaxr
##' @import stats
##' @import utils
##' @importFrom SummarizedExperiment assays colData SummarizedExperiment
##' @importFrom stringr str_order
##' @importFrom S4Vectors DataFrame
##' @importFrom DescTools DoCall
##' @importFrom MatrixExtra tcrossprod crossprod rbind_csr as.csc.matrix
##' @export
##' @md


IFAA <- function(experiment_dat,
                 testCov = NULL,
                 ctrlCov = NULL,
                 sampleIDname = NULL,
                 testMany = TRUE,
                 ctrlMany = FALSE,
                 nRef = 40,
                 nRefMaxForEsti = 2,
                 refTaxa = NULL,
                 adjust_method = "BY",
                 fdrRate = 0.05,
                 paraJobs = NULL,
                 bootB = 500,
                 standardize = FALSE,
                 sequentialRun = FALSE,
                 refReadsThresh = 0.2,
                 taxDropThresh = 0,
                 SDThresh = 0.05,
                 SDquantilThresh = 0,
                 balanceCut = 0.2,
                 verbose = TRUE) {
  allFunc <- allUserFunc()

  results <- list()
  start.time <- proc.time()[3]
  assay_name <- names(SummarizedExperiment::assays(experiment_dat))
  MicrobData <- data.frame(t(SummarizedExperiment::assays(experiment_dat)[[assay_name]]))

  MicrobData$ID_char_1234 <- rownames(MicrobData)
  CovData <- data.frame(SummarizedExperiment::colData(experiment_dat))
  CovData$ID_char_1234 <- rownames(CovData)
  linkIDname <- "ID_char_1234"

  if (verbose) {
    runMeta <- metaData(
      MicrobData = MicrobData,
      CovData = CovData,
      linkIDname = linkIDname,
      sampleIDname = sampleIDname,
      testCov = testCov,
      ctrlCov = ctrlCov,
      testMany = testMany,
      ctrlMany = ctrlMany,
      taxDropThresh = taxDropThresh,
      standardize = standardize
    )
  } else {
    runMeta <- suppressMessages(
      metaData(
        MicrobData = MicrobData,
        CovData = CovData,
        linkIDname = linkIDname,
        sampleIDname = sampleIDname,
        testCov = testCov,
        ctrlCov = ctrlCov,
        testMany = testMany,
        ctrlMany = ctrlMany,
        taxDropThresh = taxDropThresh,
        standardize = standardize
      )
    )
  }


  data <- runMeta$data
  covariatesData <- runMeta$covariatesData
  binaryInd <- runMeta$binaryInd
  covsPrefix <- runMeta$covsPrefix
  Mprefix <- runMeta$Mprefix
  testCovInd <- runMeta$testCovInd
  testCovInOrder <- runMeta$testCovInOrder
  testCovInNewNam <- runMeta$testCovInNewNam
  ctrlCov <- runMeta$ctrlCov
  microbName <- runMeta$microbName
  newMicrobNames <- runMeta$newMicrobNames
  results$covriateNames <- runMeta$xNames

  binaryInd_test <- testCovInd[testCovInd %in% binaryInd]


  rm(runMeta)

  if (length(refTaxa) > 0) {
    if (length(unique(refTaxa)) != length(refTaxa)) {
      message("Duplicated names in refTaxa are removed")
      refTaxa <- unique(refTaxa)
    }

    if (sum(refTaxa %in% microbName) != length(refTaxa)) {
      refTaxa <- refTaxa[refTaxa %in% microbName]
      message(
        "One or more of the specified reference taxa in phase 1 have no sequencing reads
      or are not in the data set."
      )
    }
  }


  if (nRefMaxForEsti < 2) {
    nRefMaxForEsti <- 2
    warning("Needs at least two final reference taxon for estimation.")
  }

  refTaxa_newNam <- newMicrobNames[microbName %in% refTaxa]

  if (verbose) {
    results$analysisResults <- Regulariz(
      data = data,
      testCovInd = testCovInd,
      testCovInOrder = testCovInOrder,
      testCovInNewNam = testCovInNewNam,
      microbName = microbName,
      nRef = nRef,
      nRefMaxForEsti = nRefMaxForEsti,
      binaryInd = binaryInd,
      binaryInd_test = binaryInd_test,
      covsPrefix = covsPrefix,
      Mprefix = Mprefix,
      refTaxa = refTaxa_newNam,
      paraJobs = paraJobs,
      adjust_method = adjust_method,
      fwerRate = fdrRate,
      bootB = bootB,
      sequentialRun = sequentialRun,
      allFunc = allFunc,
      refReadsThresh = refReadsThresh,
      SDThresh = SDThresh,
      SDquantilThresh = SDquantilThresh,
      balanceCut = balanceCut
    )
  } else {
    results$analysisResults <- suppressMessages(
      Regulariz(
        data = data,
        testCovInd = testCovInd,
        testCovInOrder = testCovInOrder,
        testCovInNewNam = testCovInNewNam,
        microbName = microbName,
        nRef = nRef,
        nRefMaxForEsti = nRefMaxForEsti,
        binaryInd = binaryInd,
        binaryInd_test = binaryInd_test,
        covsPrefix = covsPrefix,
        Mprefix = Mprefix,
        refTaxa = refTaxa_newNam,
        paraJobs = paraJobs,
        adjust_method = adjust_method,
        fwerRate = fdrRate,
        bootB = bootB,
        sequentialRun = sequentialRun,
        allFunc = allFunc,
        refReadsThresh = refReadsThresh,
        SDThresh = SDThresh,
        SDquantilThresh = SDquantilThresh,
        balanceCut = balanceCut
      )
    )
  }



  rm(data)





  totalTimeMins <- (proc.time()[3] - start.time) / 60
  message("The entire analysis took ", round(totalTimeMins, 2), " minutes")


  if (length(sampleIDname) > 0) {
    covariatesData <-
      merge(CovData[, c(sampleIDname, linkIDname)], covariatesData,
        by = linkIDname, all =
          FALSE
      )
  }
  covariatesData <-
    covariatesData[, !colnames(covariatesData) %in% c(linkIDname)]



  output_se_obj <-
    list(
      full_results = results$analysisResults$full_results,
      metadata = list(
        covariatesData = covariatesData,
        final_ref_taxon = results$analysisResults$fin_ref_taxon_name,
        ref_taxon_count =
          results$analysisResults$goodIndpRefTaxWithCount,
        ref_taxon_est = results$analysisResults$goodIndpRefTaxWithEst,
        totalTimeMins = totalTimeMins,
        fdrRate = fdrRate,
        adjust_method = adjust_method
      )
    )



  return(output_se_obj)
}
