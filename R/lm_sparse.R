
#
## sparse lm function
#

lm_sparse <- function(x, y,intercept=FALSE, tol = 1e-07){
  
  nobs   <- length(y)
  nvar   <- ncol(x) + intercept
  
  if(!methods::is(x,"dgCMatrix")) x <- MatrixExtra::as.csc.matrix(x)
  if (intercept) {
    x <- MatrixExtra::cbind2(1,x)
    qrX <- base::qr(x, tol = tol)
    dfr  <- nobs - qrX$rank
    keep <- qrX$pivot[seq_len(qrX$rank)]  
    rm(qrX)
    newX   <- x[ , keep] 
  }else{
    qrX1 <- Matrix::qr(x)
    qrX <- base::qr( as.matrix(Matrix::qrR(qrX1)), tol = tol )
    rm(qrX1)
    dfr  <- nobs - qrX$rank
    keep <- qrX$pivot[seq_len(qrX$rank)]  
    rm(qrX)
    newX   <- x[ , keep]
  }
  
  xnames <- colnames(x)
  rm(x)
  
  XTX    <- MatrixExtra::crossprod(newX)
  Xy     <- MatrixExtra::crossprod(newX, y)
  rm(newX)
  
  coef   <- Matrix::solve(XTX, Xy, tol = tol)
  coefficients       <- rep(NA, nvar)
  coefficients[keep] <- coef     
  
  yy      <- MatrixExtra::crossprod(y)
  rm(y)
  
  RSS     <- yy - 2 * MatrixExtra::crossprod(coef, Xy) + MatrixExtra::crossprod(coef, MatrixExtra::crossprod(XTX, coef))
  var_res <- as.numeric(RSS)/dfr
  se_coef <- rep(NA, nvar)
  inv     <- Matrix::solve( XTX, diag(nrow(XTX)), tol = tol)
  rm(XTX)
  
  se_coef[keep] <- sqrt(var_res * Matrix::diag(inv))
  t1            <- coefficients/se_coef
  p             <- 2 * pt(abs(t1), df = dfr, lower.tail = FALSE)
  rm(keep)
  
  coefMat<-data.frame(estimate  = coefficients,
                     std.error = se_coef,
                     t = t1,
                     p.value   = p)
  if(length(xnames)>0){
    if(intercept) xnames[1] <- "intercept"
    row.names(coefMat) <- xnames
  }else
  {if(intercept) row.names(coefMat) <- c("intercept",seq(nvar-1))
  }        
  coefMat
}

