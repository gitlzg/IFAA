## Overview

IFAA is a robust approach to make inference on the association of covariates with the absolute abundance (AA) of microbiome in an ecosystem. 

## Installation
```r
# or from CRAN
install.packages("IFAA", repos = "http://cran.us.r-project.org")

# install from GitHub:
devtools::install_github("gitlzg/IFAA")
```
## Usage

Use example datasets to run `IFAA()` function.
```r
# Detailed instructions on the package are provided in the vignettes and manual
library(IFAA)
library(SummarizedExperiment)

data(dataM)
data(dataC)
 
# merge microbiome and covariates data 
data_merged<-merge(dataM,dataC,by="id",all=FALSE)

# Seperate microbiome data and covariate data, drop id variable
dataM_sub<-data_merged[,colnames(dataM)[!colnames(dataM)%in%c("id")]]

dataC_sub<-data_merged[,colnames(dataC)]

# Create SummarizedExperiment object for inputs
test_dat<-SummarizedExperiment(assays=list(counts=t(dataM_sub)), colData=dataC_sub)

# run IFAA
results <- IFAA(experiment_dat = test_dat,
                 testCov = c("v1", "v2"),
                 ctrlCov = c("v3"),
                 sampleID = "id",
                 fdrRate = 0.15)
```


Once the analysis is done, you can extract the full resuts and significant regression coefficients along with 95% confidence intervals using this command:
```r
# to extract all results:
summary_res<-results$full_results

# to extract significant results:
sig_taxa=subset(summary_res,sig_ind==TRUE)
```

Use the same datasets to run `MZILN()` function.
```r
results <- MZILN(experiment_dat = test_dat,
                 refTaxa=c("rawCount11"),
                 allCov=c("v1","v2","v3"),
                 sampleID = "id",
                 fdrRate=0.15)
```
Regression results including confidence intervals for the targeted ratios can be extracted in the following way:
```r
# to extract the results for all ratios with rawCount11 as the denominator:
 summary_res<-results$full_results
 
# to extract results for the ratio of a specific taxon (e.g., rawCount45) over rawCount11:
 target_ratio=summary_res[summary_res$taxon=="rawCount45",]
 
# to extract all ratios having significant associations:
 sig_ratios=subset(summary_res,sig_ind==TRUE)
 ```

## References 
- Zhigang Li, Lu Tian, A. James O'Malley, Margaret R. Karagas, Anne G. Hoen, Brock C. Christensen, Juliette C. Madan, Quran Wu, Raad Z. Gharaibeh, Christian Jobin, Hongzhe Li (2021) IFAA: Robust association identification and Inference For Absolute Abundance in microbiome analyses. Journal of the American Statistical Association. 1595-1608

- Zhigang Li, Katherine Lee, Margaret Karagas, Juliette Madan, Anne Hoen, James O’Malley and Hongzhe Li (2018 ) Conditional regression based on a multivariate zero-inflated logistic normal model for modeling microbiome data. Statistics in Biosciences  10(3):587-608
