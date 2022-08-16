dataSparsCheck <- function(data,
                           Mprefix) {
  results <- list()

  # get the original sample size
  nSub <- nrow(data)
  MVarNamLength <- nchar(Mprefix)

  # get taxa variable names
  micros <-
    lapply(substr(colnames(data), 1, MVarNamLength), function(x) {
      grep(Mprefix, x)
    })
  microPositions <- which(micros == 1)
  rm(micros)

  taxaNames <- colnames(data)[microPositions]
  rm(microPositions)

  w <- data[, taxaNames, drop = FALSE]
  rm(data, taxaNames)
  overallSparsity <-
    round(100 * sum(w == 0) / (nrow(w) * ncol(w)), 2)
  message(
    overallSparsity,
    " percent of microbiome sequencing reads are zero"
  )

  # check zero taxa and subjects with zero taxa reads
  numTaxaNoReads <- length(which(Matrix::colSums(w) == 0))
  if (numTaxaNoReads > 0) {
    message(
      "There are ",
      numTaxaNoReads,
      " taxa without any sequencing reads and
        excluded from the analysis"
    )
  }
  rm(numTaxaNoReads)

  numSubNoReads <- length(which(Matrix::rowSums(w) == 0))
  if (numSubNoReads > 0) {
    message(
      "There are ",
      numSubNoReads,
      " subjects without any sequencing reads and
        excluded from the analysis."
    )
  }
  rm(numSubNoReads)
  return(list(sum(w == 0), overallSparsity))
}
