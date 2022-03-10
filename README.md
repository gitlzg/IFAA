

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
                testCov = c("v1"),
                ctrlCov = c("v2","v3"))
```


Once the analysis is done, you can extract the significant regression coefficients along with 95% confidence intervals using this command:
```r
results$sig_results
```

The function can also take csv or tsv data files directly by reading the file directory paths using the first two arguments:
```r
M="pathToTheCsvFile/microbiomeData.csv" 
# or 
M="pathToTheTsvFile/microbiomeData.tsv"

C="pathToTheCsvFile/covariatesData.csv" 
# or 
C="pathToTheTsvFile/covariatesData.tsv"

results=IFAA(MicrobData=M,CovData=C,...)
```


Use sample datasets to run `MZILN()` function.
```r
results <- MZILN(MicrobData = dataM,
                CovData = dataC,
                 linkIDname = "id",
                 targetTaxa = "rawCount18",
                 refTaxa=c("rawCount11"),
                 allCov=c("v1","v2","v3"))
```
Regression results including confidence intervals for the targeted ratios can be extracted in the following way:
```r
results$targettaxa_result_list
```

## References 
- Zhigang Li, Lu Tian, A. James O'Malley, Margaret R. Karagas, Anne G. Hoen, Brock C. Christensen, Juliette C. Madan, Quran Wu, Raad Z. Gharaibeh, Christian Jobin, Hongzhe Li (2021) IFAA: Robust association identification and Inference For Absolute Abundance in microbiome analyses. Journal of the American Statistical Association. 1595-1608

- Zhigang Li, Katherine Lee, Margaret Karagas, Juliette Madan, Anne Hoen, James O’Malley and Hongzhe Li (2018 ) Conditional regression based on a multivariate zero-inflated logistic normal model for modeling microbiome data. Statistics in Biosciences  10(3):587-608
