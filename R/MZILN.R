##' Conditional regression for microbiome analysis based on multivariate zero-inflated logistic normal model
##'
##' @description
##' The MZILN function is for estimating and testing the associations of abundance ratios with covariates.
##' \loadmathjax
##'
##' The regression model for `MZILN()` can be expressed as follows:
##' \mjdeqn{\log\bigg(\frac{\mathcal{Y}_i^k}{\mathcal{Y}_i^{K+1}}\bigg)|\mathcal{Y}_i^k>0,\mathcal{Y}_i^{K+1}>0=\alpha^{0k}+\mathcal{X}_i^T\alpha^k+\epsilon_i^k,\hspace{0.2cm}k=1,...,K}{}
##' where
##' - \mjeqn{\mathcal{Y}_i^k}{} is the AA of taxa \mjeqn{k}{} in subject \mjeqn{i}{} in the entire
##' ecosystem.
##' - \mjeqn{\mathcal{Y}_i^{K+1}}{} is the reference taxon (specified by user).
##' - \mjeqn{\mathcal{X}_i}{} is the covariate matrix for all covariates including confounders.
##' - \mjeqn{\alpha^k}{} is the regression coefficients along with their 95% confidence intervals that will be estimated by the `MZILN()` function.
##'
##' High-dimensional \mjeqn{X_i}{} is handled by regularization.
##'
##' When using this function, most of the time, users just need to feed the first three inputs to the function: `experiment_dat`, `refTaxa` and `allCov`. All other inputs can just take their default values.
##' 
##' @param experiment_dat A SummarizedExperiment object containing microbiome data and covariates (see example on how to create a SummarizedExperiment object). The microbiome data can be
##' absolute abundance or relative abundance with each column per sample and each row per taxon/OTU/ASV (or any other unit). No imputation is needed for zero-valued data points.
##' The covariates data contains covariates and confounders with each row per sample and each column per variable. The covariates data has to be numeric or binary. Categorical variables should be converted into dummy variables.
##' @param refTaxa Denominator taxa names specified by the user for the targeted ratios. This could be a vector of names.
##' @param allCov All covariates of interest (including confounders) for estimating and testing their associations with the targeted ratios. Default is 'NULL' meaning that all covariates in covData are of interest.
##' @param sampleIDname Name of the sample ID variable in the data. In the case that the data does not have an ID variable, this can be ignored. Default is NULL.
##' @param adjust_method The adjusting method for p value adjustment. Default is "BY" for dependent FDR adjustment. It can take any adjustment method for p.adjust function in R.
##' @param fdrRate The false discovery rate for identifying taxa/OTU/ASV associated with `allCov`. Default is `0.05`.
##' @param paraJobs If `sequentialRun` is `FALSE`, this specifies the number of parallel jobs that will be registered to run the algorithm. If specified as `NULL`, it will automatically detect the cores to decide the number of parallel jobs. Default is `NULL`.
##' @param bootB Number of bootstrap samples for obtaining confidence interval of estimates for the high dimensional regression. The default is `500`.
##' @param taxDropThresh The threshold of number of non-zero sequencing reads for each taxon to be dropped from the analysis. The default is `0` which means taxon without any sequencing reads will be dropped from the analysis.
##' @param standardize This takes a logical value `TRUE` or `FALSE`. If `TRUE`, the design matrix for X will be standardized in the analyses and the results. Default is `FALSE`.
##' @param sequentialRun This takes a logical value `TRUE` or `FALSE`. Default is `TRUE`. It can be set to be "FALSE" to increase speed if there are multiple taxa in the argument 'refTaxa'.
##' @param verbose Whether the process message is printed out to the console. The default is TRUE.
##' @return A list with two elements.
##'
##' - `full_results`: The main results for MZILN containing the estimation and testing results for all associations between all taxa ratios with refTaxan being the denominator and all covariates in `allCov`. It is a dataframe with each row representing an association, and ten columns named as "ref_tax", "taxon", "cov", "estimate", "SE.est", "CI.low", "CI.up", "adj.p.value", "unadj.p.value", and "sig_ind". The columns correspond to the denominator taxon, numerator taxon, covariate name, association estimates, standard error estimates, lower bound and upper bound of the 95% confidence interval, adjusted p value, and the indicator showing whether the association is significant after multiple testing adjustment.
##'
##' - `metadata`: The metadata is a list containing total time used in minutes, FDR rate, and multiple testing adjustment method used.
##'
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
##' ## Run MZILN function
##' set.seed(123) # For full reproducibility
##'
##' results <- MZILN(experiment_dat = test_dat,
##'                 refTaxa=c("rawCount11"),
##'                 allCov=c("x1","x2","x3"),
##'                 sampleIDname="ID",
##'                 fdrRate=0.05)
##' ## to extract the results for all ratios with rawCount11 as the denominator:
##' summary_res<-results$full_results
##' ## to extract results for the ratio of a specific taxon (e.g., rawCount45) over rawCount11:
##' target_ratio=summary_res[summary_res$taxon=="rawCount45",]
##' ## to extract all of the ratios having significant associations:
##' sig_ratios=subset(summary_res,sig_ind==TRUE)


##' @references Li et al.(2018) Conditional Regression Based on a Multivariate Zero-Inflated Logistic-Normal Model for Microbiome Relative Abundance Data. Statistics in Biosciences 10(3): 587-608

##' @export
##' @md

MZILN <- function(experiment_dat,
                  refTaxa,
                  allCov = NULL,
                  sampleIDname = NULL,
                  adjust_method = "BY",
                  fdrRate = 0.05,
                  paraJobs = NULL,
                  bootB = 500,
                  taxDropThresh = 0,
                  standardize = FALSE,
                  sequentialRun = FALSE,
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
      testCov = allCov,
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
        testCov = allCov,
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
  rm(runMeta)

  binaryInd_test <- testCovInd[testCovInd %in% binaryInd]

  nRef <- length(refTaxa)
  refTaxa_newNam <- newMicrobNames[microbName %in% refTaxa]


  if (length(refTaxa) > 0) {
    if (sum(refTaxa %in% microbName) != length(refTaxa)) {
      stop(
        "One or more of the specified reference taxa have no sequencing reads
      or are not in the data set. Double check the names of the reference taxa and their
           sparsity levels."
      )
    }
  }

  if (verbose) {
    results$analysisResults <- Regulariz_MZILN(
      data = data,
      nRef = nRef,
      testCovInd = testCovInd,
      testCovInOrder = testCovInOrder,
      testCovInNewNam = testCovInNewNam,
      microbName = microbName,
      binaryInd = binaryInd,
      binaryInd_test = binaryInd_test,
      covsPrefix = covsPrefix,
      Mprefix = Mprefix,
      refTaxa = refTaxa_newNam,
      adjust_method = adjust_method,
      paraJobs = paraJobs,
      fdrRate = fdrRate,
      bootB = bootB,
      sequentialRun = sequentialRun,
      allFunc = allFunc
    )
  } else {
    results$analysisResults <- suppressMessages(
      Regulariz_MZILN(
        data = data,
        nRef = nRef,
        testCovInd = testCovInd,
        testCovInOrder = testCovInOrder,
        testCovInNewNam = testCovInNewNam,
        microbName = microbName,
        binaryInd = binaryInd,
        binaryInd_test = binaryInd_test,
        covsPrefix = covsPrefix,
        Mprefix = Mprefix,
        refTaxa = refTaxa_newNam,
        adjust_method = adjust_method,
        paraJobs = paraJobs,
        fdrRate = fdrRate,
        bootB = bootB,
        sequentialRun = sequentialRun,
        allFunc = allFunc
      )
    )
  }

  rm(data)

  if (length(sampleIDname) > 0) {
    covariatesData <-
      merge(CovData[, c(sampleIDname, linkIDname)], covariatesData, all = FALSE)
  }
  covariatesData <-
    covariatesData[, !colnames(covariatesData) %in% c(linkIDname)]

  rm(testCovInOrder, ctrlCov, microbName)

  totalTimeMins <- (proc.time()[3] - start.time) / 60
  message("The entire analysis took ", round(totalTimeMins, 2), " minutes")


  output_se_obj <-
    list(
      full_results = results$analysisResults$full_results,
      metadata = list(
        covariatesData = covariatesData,
        totalTimeMins = totalTimeMins,
        fdrRate = fdrRate,
        adjust_method = adjust_method
      )
    )

  return(output_se_obj)
}
