runBootLassoHDCI <- function(x,
                             y,
                             nfolds = 5,
                             tol = 1e-05,
                             bootB,
                             bootLassoAlpha) {
  results <- list()

  nObsAll <- length(y)
  nBeta <- ncol(x)

  qrX <- qr(x = x, tol = tol)
  keep <- qrX$pivot[seq_len(qrX$rank)]
  rm(qrX)

  # remove colinear columns
  if (length(keep) < nBeta) {
    x <- x[, keep]
  }

  beta.LPR <- rep(0, nBeta)
  betaCI.LPR <- matrix(0.2, nrow = 2, ncol = nBeta)

  bootResu <- HDCI::bootLOPR(
    x = as.matrix(x),
    y = as.vector(y),
    B = bootB,
    nfolds = nfolds,
    standardize = FALSE,
    intercept = FALSE,
    alpha = bootLassoAlpha
  )

  beta.LPR[keep] <- bootResu$Beta.LPR
  betaCI.LPR[, keep] <- bootResu$interval.LPR
  rm(bootResu, keep)

  betaCIlow.LPR <- betaCI.LPR[1, ]
  betaCIhi.LPR <- betaCI.LPR[2, ]

  results$Beta.LPR <- beta.LPR
  results$interval.LPR <- rbind(betaCIlow.LPR, betaCIhi.LPR)

  return(results)
}
