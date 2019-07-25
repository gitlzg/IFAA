# Type the following line in R to install the package IFAA
 devtools::install_github("gitlzg/IFAA")
# An example to use the package
 library(IFAA)
 
 data(dataM)
 
 data(dataC)
 
 results=IFAA(Microbdata=dataM,CovData=dataC,linkIDname="id",
 
  testCov=c("v1","v2"),ctrlCov=c("v3"),nRef=2,
      
  refTaxa=c("rawCount1"),nPermu=2,paraJobs=NULL,
      
  reguMethod="mcp",fwerRate=0.3,bootB=3)
