metaData <- function(MicrobData,
                     CovData,
                     linkIDname,
                     sampleIDname,
                     testCov = NULL,
                     ctrlCov = NULL,
                     testMany = TRUE,
                     ctrlMany = FALSE,
                     MZILN = FALSE,
                     taxDropThresh,
                     standardize) {
  results <- list()

  testCov <- unique(testCov)
  ctrlCov <- unique(ctrlCov)

  if (length(linkIDname) == 0) {
    stop("linkIDname is missing.")
  }

  if (length(testCov) > 0 | length(ctrlCov) > 0) {
    if (sum(c(testCov, ctrlCov) %in% colnames(CovData)) != length(c(testCov, ctrlCov))) {
      stop("some covariates are not available in the data.")
    }
  }

  if (sum(testCov %in% ctrlCov) > 0) {
    warning("Variables appeared in both testCov list and ctrlCov list will be treated as testCov.")
  }

  # read microbiome data
  if (is.matrix(MicrobData)) {
    MdataWithId <- data.matrix(MicrobData)
  }
  if (is.data.frame(MicrobData)) {
    MdataWithId <- data.matrix(MicrobData)
  }


  if (length(colnames(MdataWithId)) != ncol(MdataWithId)) {
    stop("Microbiome data lack variable names.")
  }

  MdataWithoutId <-
    data.matrix(MdataWithId[, !(colnames(MdataWithId) %in% c(linkIDname, sampleIDname)),
      drop =
        FALSE
    ])
  uniqMnames <- unique(colnames(MdataWithoutId))
  if (length(uniqMnames) != length(colnames(MdataWithoutId))) {
    nDup <- length(colnames(MdataWithoutId)) - length(uniqMnames)
    message(
      nDup,
      " Duplicated taxa/OTU/ASV names are removed from the microbiome data."
    )
  }
  MdataWithoutId <- MdataWithoutId[, uniqMnames]
  MdataWithId <-
    cbind(MdataWithId[, linkIDname, drop = FALSE], MdataWithoutId)

  if (!all(MdataWithoutId >= 0)) {
    stop("Microbiome data contains negative values.")
  }
  rm(MdataWithoutId)

  missPropMData <-
    sum(is.na(MdataWithId[, linkIDname])) / nrow(MdataWithId)
  if (missPropMData > 0.8) {
    warning(
      "There are over 80% missing values for the linkId variable in the Microbiome data file.
             Double check the data format."
    )
  }

  # read covariate data
  if (is.matrix(CovData)) {
    CovarWithId <- data.matrix(CovData)
  }
  if (is.data.frame(CovData)) {
    CovarWithId <- data.matrix(CovData)
  }


  if (length(colnames(CovarWithId)) != ncol(CovarWithId)) {
    stop("Covariate data lack variable names.")
  }

  missPropCovData <-
    sum(is.na(CovarWithId[, linkIDname])) / nrow(CovarWithId)
  if (missPropCovData > 0.8) {
    warning(
      "There are over 80% missing values for the linkId variable in the covariates data file.
             Double check the data format."
    )
  }

  Covariates1 <-
    CovarWithId[, !colnames(CovarWithId) %in% c(linkIDname, sampleIDname),
      drop =
        FALSE
    ]

  # determine testCov and ctrlCov
  if (length(testCov) == 0) {
    if (!testMany) {
      stop("No covariates are specified for estimating associations of interest.")
    } else {
      message(
        "Associations are being estimated for all covariates since no covariates are specified for testCov."
      )
      testCov <- colnames(Covariates1)
    }
  }
  results$testCov <- testCov

  xNames <- colnames(Covariates1)
  rm(Covariates1)

  if (length(ctrlCov) == 0 & ctrlMany) {
    message("All variables except testCov are considered as control covariates.")
    ctrlCov <- xNames[!xNames %in% testCov]
  }

  # make sure testCov and ctrlCov mutually exclusive
  ctrlCov <- ctrlCov[!(ctrlCov %in% testCov)]
  results$ctrlCov <- ctrlCov

  # merge data to remove missing
  CovarWithId1 <- CovarWithId[, c(linkIDname, testCov, ctrlCov)]

  allRawData <- data.matrix(na.omit(
    merge(
      CovarWithId1,
      MdataWithId,
      by = linkIDname,
      all.x = FALSE,
      all.y = FALSE
    )
  ))

  CovarWithId <-
    allRawData[, (colnames(allRawData) %in% colnames(CovarWithId1)),
      drop =
        FALSE
    ]
  Covariates <-
    CovarWithId[, !colnames(CovarWithId) %in% linkIDname,
      drop =
        FALSE
    ]
  rm(CovarWithId1)


  if (!is.numeric(Covariates[, testCov, drop = FALSE])) {
    warning("There are non-numeric variables in the test covariates")
    nTestCov <- length(testCov)
    numCheck <- unlist(lapply(seq(nTestCov), function(i) {
      is.numeric(Covariates[, testCov[i]])
    })) + 0
    for (i in which(numCheck == 0)) {
      Covariates[, testCov[i]] <-
        as.numeric(factor(Covariates[, testCov[i]]))
    }
  }


  MdataWithId <-
    allRawData[, (colnames(allRawData) %in% colnames(MdataWithId))]
  Mdata_raw <-
    MdataWithId[, !(colnames(MdataWithId) %in% linkIDname),
      drop =
        FALSE
    ]
  rm(allRawData)

  # check zero taxa and subjects with zero taxa reads
  numTaxaNoReads <- sum(colSums(Mdata_raw != 0) <= taxDropThresh)
  if (numTaxaNoReads > 0) {
    Mdata_raw <- Mdata_raw[, !(colSums(Mdata_raw != 0) <= taxDropThresh)]
    message(
      "There are ",
      numTaxaNoReads,
      " taxa without any sequencing reads before
        data merging, and excluded from the analysis"
    )
  }
  rm(numTaxaNoReads)

  numSubNoReads <- sum(rowSums(Mdata_raw != 0) <= 1)
  if (numSubNoReads > 0) {
    message(
      "There are ",
      numSubNoReads,
      " samples with zero or one sequencing read and
        excluded from the analysis"
    )
    subKeep <- !(rowSums(Mdata_raw != 0) <= 1)
    Mdata_raw <- Mdata_raw[subKeep, ]
    MdataWithId <- MdataWithId[subKeep, ]
    rm(subKeep)
  }
  rm(numSubNoReads)

  Mdata <- Mdata_raw
  rm(Mdata_raw)

  microbName1 <- colnames(Mdata)
  microbName <- microbName1
  newMicrobNames1 <- paste0("microb", seq(length(microbName)))
  newMicrobNames <- newMicrobNames1
  results$Mprefix <- "microb"

  colnames(Mdata) <- newMicrobNames
  MdataWithId_new <-
    cbind(MdataWithId[, linkIDname, drop = FALSE], Mdata)
  results$microbName <- microbName
  results$newMicrobNames <- newMicrobNames

  if (sum(is.na(Covariates)) > 0) {
    message("Samples with missing covariate values are removed from the analysis.")
  }

  if (!is.numeric(Covariates[, ctrlCov, drop = FALSE])) {
    warning("There are non-numeric variables in the control covariates")
    nCtrlCov <- length(ctrlCov)
    numCheck <- unlist(lapply(seq(nCtrlCov), function(i) {
      is.numeric(Covariates[, ctrlCov[i]])
    })) + 0
    for (i in which(numCheck == 0)) {
      Covariates[, ctrlCov[i]] <-
        as.numeric(factor(Covariates[, ctrlCov[i]]))
    }
  }


  xNames <- colnames(Covariates)
  nCov <- length(xNames)

  binCheck <- unlist(lapply(seq(nCov), function(i) {
    dim(table(Covariates[, xNames[i]]))
  }))

  if ((sum(binCheck == 2)) > 0) {
    binaryInd <- which(binCheck == 2)
    results$varNamForBin <- xNames[binCheck == 2]
    results$nBinVars <- length(results$varNamForBin)
    for (i in results$varNamForBin) {
      mini <- min(Covariates[, i], na.rm = TRUE)
      maxi <- max(Covariates[, i], na.rm = TRUE)
      if (!(mini == 0 & maxi == 1)) {
        Covariates[Covariates[, i] == mini, i] <- 0
        Covariates[Covariates[, i] == maxi, i] <- 1
        message(
          "Binary covariate",
          i,
          "is not coded as 0/1 which may generate
                  analysis bias. It has been changed to 0/1. The changed covariates
                  data can be extracted from the result file."
        )
      }
    }
  } else {
    results$nBinVars <- 0
    binaryInd <- NULL
    results$varNamForBin <- NULL
  }

  results$binaryInd <- binaryInd
  results$xNames <- colnames(Covariates)
  if (standardize) {
    if (length(binaryInd) > 0) {
      Covariates[, -binaryInd] <- Covariates[, -binaryInd] %*%
        diag(1 / apply(Covariates[, -binaryInd], 2, function(x) {
          sd(x, na.rm = TRUE)
        }))
    } else {
      Covariates <- Covariates %*%
        diag(1 / apply(Covariates, 2, function(x) {
          sd(x, na.rm = TRUE)
        }))
    }
  }

  #### Test for constant column
  sd_col <- apply(Covariates, 2, function(x) {
    sd(x, na.rm = TRUE)
  })
  sd_zero_loc <- which(sd_col == 0)
  if (length(sd_zero_loc) > 0) {
    Covariates[, sd_zero_loc] <- 0
    warning(
      "Covariate ",
      colnames(Covariates)[sd_zero_loc],
      " is constant, please double check"
    )
  }

  xNewNames <- paste0("x", seq(length(xNames)))
  colnames(Covariates) <- xNewNames
  results$covsPrefix <- "x"
  results$xNewNames <- xNewNames

  results$testCovInd <- which((results$xNames) %in% testCov)
  results$testCovInOrder <- results$xNames[results$testCovInd]
  results$testCovInNewNam <- results$xNewNames[results$testCovInd]
  rm(xNames, xNewNames)

  CovarWithId_new <-
    cbind(CovarWithId[, linkIDname, drop = FALSE], Covariates)

  data <- merge(
    MdataWithId_new,
    CovarWithId_new,
    by = linkIDname,
    all.x = FALSE,
    all.y = FALSE
  )
  dataOmit <- na.omit(data)

  results$covariatesData <- CovarWithId_new
  colnames(results$covariatesData) <- c(linkIDname, results$xNames)
  rm(MdataWithId_new, CovarWithId_new)

  Mdata_omit <- dataOmit[, newMicrobNames1]

  # check taxa with zero or 1 read again after all missing data removed
  numTaxaNoReads <- sum(colSums(Mdata_omit != 0) <= taxDropThresh)
  if (numTaxaNoReads == 0) {
    results$data <- data.matrix(dataOmit)
  }
  if (numTaxaNoReads > 0) {
    dataOmit_noTaxa <-
      dataOmit[, !(colnames(dataOmit) %in% newMicrobNames1)]
    microbToRetain <-
      newMicrobNames1[!(colSums(Mdata_omit != 0) <= taxDropThresh)]
    message(
      "There are ",
      numTaxaNoReads,
      " taxa without any sequencing reads after merging
        and removing all missing data, and excluded from the analysis"
    )
    MdataToRetain <- Mdata_omit[, microbToRetain]
    microbName <- microbName1[newMicrobNames1 %in% microbToRetain]
    results$microbName <- microbName
    newMicrobNames1 <- paste0("microb", seq(length(microbName)))
    newMicrobNames <- newMicrobNames1
    results$newMicrobNames <- newMicrobNames
    colnames(MdataToRetain) <- microbToRetain
    results$data <-
      data.matrix(cbind(dataOmit_noTaxa, MdataToRetain))
    rm(dataOmit_noTaxa, MdataToRetain, microbToRetain)
  }
  rm(
    numTaxaNoReads,
    data,
    microbName,
    microbName1,
    newMicrobNames,
    newMicrobNames1
  )

  # output data summary
  message("Data dimensions (after removing missing data if any):")
  message(dim(results$data)[1], " samples")
  message(ncol(Mdata), " taxa/OTU/ASV")

  if (!MZILN) {
    message(
      length(results$testCovInOrder),
      " testCov variables in the analysis"
    )
  }
  if (MZILN) {
    message(length(results$testCovInOrder), " covariates in the analysis")
  }

  if (length(results$testCovInOrder) > 0) {
    if (!MZILN) {
      message("These are the testCov variables:")
    }
    if (MZILN) {
      message("These are the covariates:")
    }
    testCovPrint <- results$testCovInOrder[1]
    if (length(results$testCovInOrder) > 1) {
      for (i in 2:length(results$testCovInOrder)) {
        testCovPrint <-
          paste0(testCovPrint, ", ", results$testCovInOrder[i])
      }
    }
    message(testCovPrint)
  }
  rm(testCov)
  if (!MZILN) {
    message(length(ctrlCov), " ctrlCov variables in the analysis ")
    if (length(ctrlCov) > 0) {
      message("These are the ctrlCov variables:")
      ctrlCovPrint <- ctrlCov[1]
      if (length(ctrlCov) > 1) {
        for (i in 2:length(ctrlCov)) {
          ctrlCovPrint <- paste0(ctrlCovPrint, ", ", ctrlCov[i])
        }
      }
      message(ctrlCovPrint)
    }
    rm(ctrlCov)
  }
  message(results$nBinVars, " binary covariates in the analysis")
  if (results$nBinVars > 0) {
    message("These are the binary covariates:")
    out_message <- paste0(results$varNamForBin, " ")
    message(out_message)
  }
  rm(Mdata, Covariates, binCheck)
  return(results)
}
