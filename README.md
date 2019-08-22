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
      
  It takes about 21 minutes to finish the example analysis on a 8-core machine.
  
  Once the analysis is done, you extract the regression coefficients along with 95% confidence intervals using this command:
  
  results$analysisResults$estByCovList
  
