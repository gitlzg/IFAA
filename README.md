## Overview

IFAA offers a robust approach to make inference on the associations of covariates 
with the absolute abundance (AA) of microbiome in an ecosystem, and the associations of covariates with the abundance ratios of microbiome taxa/OTU/ASV (or other units). 

## Installation
```r
# Install from CRAN
install.packages("IFAA", repos = "http://cran.us.r-project.org")

# Or from GitHub:
devtools::install_github("gitlzg/IFAA")
```
## Usage

This is a brief example to run `IFAA()` and `MZILN()` functions. More detailed instructions can be found in the manual and vignettes files. This example generates a toy  example data from scratch, with 10 taxon, 20 subjects, and 3 covariates.
```r
# Detailed instructions on the package are 
# provided in the vignettes and manual
library(IFAA)

## create an ID variable for the example data
ID=seq_len(20)

## generate three covariates x1, x2, and x3, with x2 binary
x1<-rnorm(20)
x2<-rbinom(20,1,0.5)
x3<-rnorm(20)
dataC<-data.frame(cbind(ID,x1,x2,x3))

## Coefficients for x1, x2, and x3 among 10 taxa.
beta_1<-c(0.1,rep(0,9))
beta_2<-c(0,0.2,rep(0,8))
beta_3<-rnorm(10)
beta_mat<-cbind(beta_1,beta_2,beta_3)

## Generate absolute abundance for 10 taxa in ecosystem.
dataM_eco<-floor(exp(10+as.matrix(dataC[,-1]) %*% t(beta_mat) + rnorm(200,sd=0.05)))

## Generate sequence depth and generate observed abundance
Ci<-runif(20,0.01,0.05)
dataM<-floor(apply(dataM_eco,2,function(x) x*Ci))
colnames(dataM)<-paste0("rawCount",1:10)

## Randomly introduce 0 to make 25% sparsity level.
dataM[sample(seq_len(length(dataM)),length(dataM)/4)]<-0

dataM<-data.frame(cbind(ID,dataM))

## The following steps are to create a SummarizedExperiment object.
## If you already have a SummarizedExperiment format data, you can
## ignore the following steps and directly feed it to the IFAA function.

## Merge two dataset by ID variable
data_merged<-merge(dataM,dataC,by="ID",all=FALSE)

## Seperate microbiome data and covariate data, drop ID variable from microbiome data
dataM_sub<-data_merged[,colnames(dataM)[!colnames(dataM)%in%c("ID")]]
dataC_sub<-data_merged[,colnames(dataC)]

## Create SummarizedExperiment object
test_dat<-SummarizedExperiment::SummarizedExperiment(
assays=list(MicrobData=t(dataM_sub)), colData=dataC_sub)

# Again, if you already have a SummarizedExperiment format data, you can 
# ignore the above steps

# run IFAA
set.seed(123) # For full reproducibility
results <- IFAA(experiment_dat = test_dat,
                testCov = c("x1","x2"),
                ctrlCov = c("x3"),
                sampleIDname="ID",
                fdrRate = 0.05)
```


Once the analysis is done, you can extract the full resuts and significant regression coefficients along with 95% confidence intervals using this command:
```r
# to extract all results:
summary_res<-results$full_results

# to extract significant results:
sig_results=subset(summary_res,sig_ind==TRUE)
```

Use the same datasets to run `MZILN()` function.
```r
set.seed(123) # For full reproducibility
resultsRatio <- MZILN(
  experiment_dat = test_dat,
  microbVar = c("rawCount5","rawCount8"),
  refTaxa = c("rawCount10"),
  allCov = c("x1", "x2", "x3"),
  sampleIDname = "ID",
  fdrRate = 0.05,
  paraJobs = 2
)
```
Regression results including confidence intervals for the targeted ratios can be extracted in the following way:
```r
# to extract results for the ratio "rawCount5/rawCount10":
 target_ratio5=summary_res[summary_res$taxon=="rawCount5",]

# to extract results for the ratio "rawCount8/rawCount10":
 target_ratio8=summary_res[summary_res$taxon=="rawCount8",]
 ```
 

## References 
- Zhigang Li, Lu Tian, A. James O'Malley, Margaret R. Karagas, Anne G. Hoen, Brock C. Christensen, Juliette C. Madan, Quran Wu, Raad Z. Gharaibeh, Christian Jobin, Hongzhe Li (2021) IFAA: Robust association identification and Inference For Absolute Abundance in microbiome analyses. Journal of the American Statistical Association. 1595-1608

- Zhigang Li, Katherine Lee, Margaret Karagas, Juliette Madan, Anne Hoen, James O'Malley and Hongzhe Li (2018) Conditional regression based on a multivariate zero-inflated logistic normal model for modeling microbiome data. Statistics in Biosciences 10(3):587-608