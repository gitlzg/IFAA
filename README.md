# IFAA 

## Overview

IFAA is a novel approach to make inference on the association of covariates with the absolute abundance (AA) of microbiome in an ecosystem. 

## Installation
```r
# install from GitHub:
devtools::install_github("gitlzg/IFAA")
```
## Usage

Use sample datasets to run `IFAA()` function.
```r
# Detailed instructions on the package are provided in the manual and vignette
library(IFAA)
 
data(dataM)
data(dataC)
 
results <- IFAA(MicrobData = dataM,
                CovData = dataC,
                linkIDname = "id",
                testCov = c("v1", "v2"),
                ctrlCov = c("v3"),
                nRef = 4,
                nPermu = 4,
                bootB = 5)
```


Once the analysis is done, you can extract the regression coefficients along with 95% confidence intervals using this command:
```r
results$analysisResults$estByCovList
```

The function can also take csv or tsv data files directly by reading the file directory paths using the first two arguments:
```r
M="pathToTheCsvFile//microbiomeData.csv" 
# or 
M="pathToTheTsvFile//microbiomeData.tsv"

C="pathToTheCsvFile//covariatesData.csv" 
# or 
C="pathToTheTsvFile//covariatesData.tsv"

results=IFAA(MicrobData=M,CovData=C,...)
```


Use sample datasets to run `MZILN()` function.
```r
results <- MZILN(MicrobData = dataM,
                 CovData = dataC,
                 linkIDname = "id",
                 allCov = c("v1","v2","v3"),
                 refTaxa=c("rawCount11")
                 )
```
Regression results including confidence intervals can be extracted in the following way:
```r
results$analysisResults$estByRefTaxaList$rawCount11$estByCovList
```

## References 
- Li Z, Tian L, O'Malley J, Karagas MR, Hoen AG, Christensen BC, Madan JC, Wu Q, Gharaibeh RZ, Jobin C, Li H (In press) IFAA: Robust association identification and Inference For Absolute Abundance in microbiome analyses. Journal of the American Statistical Association

- Li Z, Lee K, Karagas MR, Madan JC, Hoen AG, O'Malley J and Li H (2018 ) Conditional regression based on a multivariate zero-inflated logistic normal model for modeling microbiome data. Statistics in Biosciences  10(3):587-608
