dataRecovTrans <- function(data,
                           ref,
                           Mprefix,
                           covsPrefix,
                           binPredInd,
                           contCovStd = FALSE,
                           xOnly = FALSE,
                           yOnly = FALSE) {
  results <- list()

  # load A and log-ratio transformed RA
  data.and.init <- AIcalcu(
    data = data,
    ref = ref,
    Mprefix = Mprefix,
    covsPrefix = covsPrefix,
    binPredInd = binPredInd,
    contCovStd = contCovStd
  )

  rm(data)

  taxaNames <- data.and.init$taxaNames
  A <- data.and.init$A
  logRatiow <- data.and.init$logRatiow
  nSub <- data.and.init$nSub
  nTaxa <- data.and.init$nTaxa
  xData <- data.and.init$xData
  nPredics <- data.and.init$nPredics
  twoList <- data.and.init$twoList
  lLast <- data.and.init$lLast
  L <- data.and.init$l
  lengthTwoList <- data.and.init$lengthTwoList
  rm(data.and.init)

  nNorm <- nTaxa - 1
  xDimension <- nPredics + 1 # predictors+intercept

  # create omegaRoot
  omegaRoot <- list()
  for (j in seq_len(lengthTwoList)) {
    i <- twoList[j]
    if (lLast[i] == nTaxa) {
      omegaRoot[[i]] <- Matrix::Diagonal(L[i] - 1)
    } else {
      if (L[i] == 2) {
        omegaRoot[[i]] <- sqrt(0.5)
      } else {
        dim <- L[i] - 1
        a <- (1 + (dim - 2) / 2) / (1 / 2 * (1 + (dim - 1) / 2))
        b <- -1 / (1 + (dim - 1) / 2)

        # calculate the square root of omega assuming it is exchangeable
        aStar <- dim^2 / ((dim - 1)^2)
        bStar <-
          b * (dim - 2) / (dim - 1) - a * (dim^2 - 2 * dim + 2) /
            ((dim - 1)^2)
        cStar <- (0.5 * b - 0.5 * a * (dim - 2) / (dim - 1))^2
        cSquare <-
          (-bStar + sqrt(bStar^2 - 4 * aStar * cStar)) / (2 * aStar)
        if (cSquare < 0) {
          stop("no solution for square root of omega")
        }
        d <- sqrt((0.5 * a - cSquare) / (dim - 1))
        if (is.na(d)) {
          stop("no solution for off-diagnal elements for square root of omega")
        }
        c <- (0.5 * b - (dim - 2) * (d^2)) / (2 * d)
        omegaRoot[[i]] <-
          -((c - d) * Matrix::Diagonal(dim) + d * matrix(1,
            nrow =
              dim, ncol = dim
          ))
      }
    }
  }
  rm(L, lLast)

  if (xOnly) {
    # create X_i in the regression equation using Kronecker product
    xDataWithInter <- as.matrix(cbind(rep(1, nSub), xData))
    rm(xData)

    xtild <- list()
    for (j in seq_len(lengthTwoList)) {
      i <- twoList[j]
      non0v <- rep(xDataWithInter[i, ], nNorm)
      posi1 <- seq(nPredics + 1)
      posi2 <- matrix(data = (rep((
        seq(0, (nNorm - 1)) * (nPredics + 1) * nNorm
      ), each = (nPredics + 1))), nrow = (nPredics + 1))
      posi3 <- matrix(data = (rep((
        seq(0, (nNorm - 1)) * (nPredics + 1)
      ), each = (nPredics + 1))), nrow = (nPredics + 1))
      posi <- c(posi1 + posi2 + posi3)
      rm(posi2, posi3)

      sVec <- Matrix::sparseVector(
        x = non0v,
        i = posi,
        length = ((nPredics + 1) * nNorm * nNorm)
      )

      xInRegres.i.transpos <-
        Matrix::Matrix(sVec,
          nrow = (nNorm * (nPredics + 1)), ncol =
            nNorm
        )
      # rm(sVec)

      xDataTilda.i <-
        MatrixExtra::crossprod(omegaRoot[[i]], MatrixExtra::tcrossprod(A[[i]], xInRegres.i.transpos))
      rm(xInRegres.i.transpos)
      xtild[[j]] <- xDataTilda.i

      rm(xDataTilda.i)
    }

    xTildalong <- DescTools::DoCall(MatrixExtra::rbind_csr, xtild)

    rm(xtild, xDataWithInter)

    results$xTildalong <- xTildalong
    rm(xTildalong)
    return(results)
  }

  if (yOnly) {
    for (j in seq_len(lengthTwoList)) {
      i <- twoList[j]
      Utilda.i <- MatrixExtra::crossprod(omegaRoot[[i]], logRatiow[[i]])
      if (j == 1) {
        UtildaLong <- as.numeric(Utilda.i)
      } else {
        UtildaLong <- c(UtildaLong, as.numeric(Utilda.i))
      }
    }
    rm(omegaRoot, logRatiow, Utilda.i)

    results$UtildaLong <- UtildaLong
    rm(UtildaLong)
    return(results)
  }

  # create X_i in the regression equation using Kronecker product
  xDataWithInter <- data.matrix(cbind(rep(1, nSub), xData))
  colnames(xDataWithInter)[1] <- "Inter"
  rm(xData)

  xtild <- list()
  for (j in seq_len(lengthTwoList)) {
    i <- twoList[j]
    non0v <- rep(xDataWithInter[i, ], nNorm)
    posi1 <- seq(nPredics + 1)
    posi2 <- matrix(data = (rep((
      seq(0, (nNorm - 1)) * (nPredics + 1) * nNorm
    ), each = (nPredics + 1))), nrow = (nPredics + 1))
    posi3 <- matrix(data = (rep((
      seq(0, (nNorm - 1)) * (nPredics + 1)
    ), each = (nPredics + 1))), nrow = (nPredics + 1))
    posi <- c(posi1 + posi2 + posi3)
    rm(posi2, posi3)

    sVec <- Matrix::sparseVector(
      x = non0v,
      i = posi,
      length = ((nPredics + 1) * nNorm * nNorm)
    )

    xInRegres.i.transpos <-
      Matrix::Matrix(sVec,
        nrow = (nNorm * (nPredics + 1)), ncol =
          nNorm
      )
    # rm(sVec)

    xDataTilda.i <-
      MatrixExtra::crossprod(omegaRoot[[i]], MatrixExtra::tcrossprod(A[[i]], xInRegres.i.transpos))
    rm(xInRegres.i.transpos)
    xtild[[j]] <- xDataTilda.i

    rm(xDataTilda.i)
  }

  xTildalong <- DescTools::DoCall(MatrixExtra::rbind_csr, xtild)

  rm(xtild, xDataWithInter)

  for (j in seq_len(lengthTwoList)) {
    i <- twoList[j]
    Utilda.i <- MatrixExtra::crossprod(omegaRoot[[i]], logRatiow[[i]])
    if (j == 1) {
      UtildaLong <- as.numeric(Utilda.i)
    } else {
      UtildaLong <- c(UtildaLong, as.numeric(Utilda.i))
    }
  }

  rm(omegaRoot, logRatiow, Utilda.i)

  # return objects
  results$UtildaLong <- UtildaLong
  rm(UtildaLong)
  results$xTildalong <- xTildalong
  rm(xTildalong)
  results$taxaNames <- taxaNames
  rm(taxaNames)
  return(results)
}
