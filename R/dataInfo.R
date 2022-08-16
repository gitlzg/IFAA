dataInfo <- function(data,
                     Mprefix,
                     covsPrefix,
                     qualifyRefTax = FALSE,
                     binPredInd,
                     refReadsThresh,
                     SDThresh,
                     SDquantilThresh,
                     balanceCut) {
  results <- list()

  # get the original sample size
  nSub <- nrow(data)
  MVarNamLength <- nchar(Mprefix)
  # get taxa variable names
  micros <- lapply(substr(colnames(data), 1, MVarNamLength), function(x) {
    grep(Mprefix, x)
  })
  microPositions <- which(micros == 1)
  rm(micros)

  nTaxa <- length(microPositions)
  taxaNames <- colnames(data)[microPositions]
  rm(microPositions)

  if (qualifyRefTax) {
    qualifyData <- data[rowSums(data[, taxaNames] > 0) >= 2, , drop = FALSE]
    w <- qualifyData[, taxaNames]
    nSubQualif <- nrow(qualifyData)
    taxaOverThresh <- taxaNames[(Matrix::colSums(w > 0) >= (nSubQualif *
      refReadsThresh))]

    # check the sd threshold
    sdTaxaOverThresh <- rep(0, length(taxaOverThresh))
    for (i in seq_len(length(taxaOverThresh))) {
      taxa.i <- w[, taxaOverThresh[i]]
      if (sum(taxa.i > 0) > 1) {
        sdTaxaOverThresh[i] <- sd(taxa.i[(taxa.i > 0)])
      }
    }

    results$sdTaxa <- sdTaxaOverThresh

    TaxaOverSdThresh <- taxaOverThresh[(sdTaxaOverThresh >= SDThresh)]
    sdOverSdThresh <- sdTaxaOverThresh[sdTaxaOverThresh >= SDThresh]

    rm(taxa.i, taxaOverThresh)

    # check the sd quantile threshold
    sdAllTaxa <- rep(0, nTaxa)
    for (i in seq_len(nTaxa)) {
      taxaAll.i <- w[, taxaNames[i]]
      posTaxaAll.i <- taxaAll.i[(taxaAll.i > 0)]
      if (length(posTaxaAll.i) > 1) {
        sdAllTaxa[i] <- sd(posTaxaAll.i)
      }
    }
    goodRefTaxaCandi <- TaxaOverSdThresh[(sdOverSdThresh >= quantile(sdAllTaxa,
      probs =
        SDquantilThresh
    ))]
    rm(sdAllTaxa, posTaxaAll.i, TaxaOverSdThresh)

    rm(w)
  }

  # get predictor data
  xVarNamLength <- nchar(covsPrefix)

  predics <- lapply(substr(colnames(data), 1, xVarNamLength), function(x) {
    grep(covsPrefix, x)
  })

  predPositions <- which(predics == 1)
  predNames <- colnames(data)[predPositions]
  nPredics <- length(predNames)
  rm(predics, predPositions, data)

  if (qualifyRefTax) {
    # find the pairs of binary preds and taxa for which the assocaiton is not identifiable
    if (length(binPredInd) > 0) {
      allBinPred <- predNames[binPredInd]
      nBinPred <- length(allBinPred)

      taxaBalanceBin <- c()

      for (i in seq_len(nTaxa)) {
        for (j in seq_len(nBinPred)) {
          twoColumns.ij <- qualifyData[, c(taxaNames[i], allBinPred[j])]
          nNonZero <- length(which(twoColumns.ij[, 1] > 0))
          sumOfBin <- sum(twoColumns.ij[(twoColumns.ij[, 1] > 0), 2])
          rm(twoColumns.ij)
          if (min(sumOfBin, (nNonZero - sumOfBin)) >= floor(balanceCut *
            nSubQualif)) {
            taxaBalanceBin <- c(taxaBalanceBin, taxaNames[i])
          }
        }
      }

      rm(allBinPred, qualifyData)
      taxaBalanceBin <- unique(taxaBalanceBin)

      # keep balanced taxa
      goodRefTaxaCandi <- goodRefTaxaCandi[(goodRefTaxaCandi %in% taxaBalanceBin)]
    }
    results$goodRefTaxaCandi <- goodRefTaxaCandi
    rm(goodRefTaxaCandi)
  }
  # return
  results$taxaNames <- taxaNames
  rm(taxaNames)
  results$predNames <- predNames
  rm(predNames)
  results$nTaxa <- nTaxa
  results$nSub <- nSub
  results$nPredics <- nPredics
  return(results)
}
