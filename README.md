# Type the following line in R to install the package IFAA
 devtools::install_github("gitlzg/IFAA")

# An example to use the package
library(IFAA)
 
data(dataM)
 
data(dataC)
 
results=IFAA(MicrobData=dataM,
CovData=dataC,
linkIDname="id",
testCov=c("v1","v2"),
ctrlCov=c("v3"),
nRef=4,
nPermu=4,
bootB=5)
      
# Once the analysis is done, you can extract the regression coefficients along with 95% confidence intervals using this command:
results$analysisResults$estByCovList
  
# The function can also take csv or tsv data files directly by reading the file directory paths using the first two arugments:
M="pathToTheCsvFile/microbiomeData.csv" or M="pathToTheTsvFile/microbiomeData.tsv"

C="pathToTheCsvFile/covariatesData.csv" or C="pathToTheTsvFile/covariatesData.tsv"

results=IFAA(MicrobData=M,CovData=C,...)

# Reference: 
Li et al. IFAA: Robust association identification and Inference For Absolute Abundance in microbiome analyses. https://arxiv.org/abs/1909.10101
