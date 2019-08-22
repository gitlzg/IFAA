# Type the following line in R to install the package IFAA
 devtools::install_github("gitlzg/IFAA")

# An example to use the package
library(IFAA)
 
data(dataM)
 
data(dataC)
 
results=IFAA(Microbdata=dataM,
 
CovData=dataC,

linkIDname="id",

testCov=c("v1","v2"),

ctrlCov=c("v3"),

nRef=40,

nPermu=50,

reguMethod="mcp",

fwerRate=0.2,

bootB=500)
      
# It takes about 21 minutes to finish the example analysis on a 8-core machine.
  
# Once the analysis is done, you can extract the regression coefficients along with 95% confidence intervals using this command:
results$analysisResults$estByCovList
  
# The function can also directly read csv or tsv files:
M="pathToTheCsvFile/microbiomeData.csv" or M="pathToTheTsvFile/microbiomeData.tsv"

C="pathToTheCsvFile/covariatesData.csv" or C="pathToTheTsvFile/covariatesData.tsv"

results=IFAA(Microbdata=M,
 
CovData=C,

linkIDname="id",

testCov=c("v1","v2"),

ctrlCov=c("v3"),

nRef=40,

nPermu=50,

reguMethod="mcp",

fwerRate=0.2,

bootB=500)
      
