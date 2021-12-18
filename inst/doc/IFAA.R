## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(IFAA)

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("IFAA", repos = "http://cran.us.r-project.org")

## ----eval=FALSE---------------------------------------------------------------
#  require(devtools)
#  devtools::install_github("gitlzg/IFAA")

## -----------------------------------------------------------------------------
library(IFAA)

data(dataM)
dim(dataM)
dataM[1:5, 1:8]

data(dataC)
dim(dataC)
dataC[1:5, ]

## ---- eval=T------------------------------------------------------------------
results <- IFAA(MicrobData = dataM,
                CovData = dataC,
                linkIDname = "id",
                testCov = c("v1", "v2"),
                ctrlCov = c("v3"),
                nRef = 3,
                paraJobs = 2,
                fdrRate = 0.25)

## ----eval=T-------------------------------------------------------------------
results$analysisResults$sig_list_each_mean

## ----eval=T-------------------------------------------------------------------
results$covariatesData

## -----------------------------------------------------------------------------
data(dataM)
dim(dataM)
dataM[1:5, 1:8]

data(dataC)
dim(dataC)
dataC[1:5, ]

## ---- eval=T------------------------------------------------------------------
results <- MZILN(MicrobData = dataM,
                CovData = dataC,
                 linkIDname = "id",
                 allCov=c("v1","v2","v3"),
                 targetTaxa = "rawCount6",
                 refTaxa=c("rawCount11"),
                 paraJobs=2)

## ----eval=T-------------------------------------------------------------------
results$analysisResults$targettaxa_result_list

## ----eval=T-------------------------------------------------------------------
results$covariatesData

