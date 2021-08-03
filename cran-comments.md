## Test environments
* local  Windows 10 x64 (build 18362), release R 4.0.2, 
* Ubuntu 16.04.6 LTS (on travis-ci), release R 4.0.2, oldrel R 3.6.3, devel (2020-10-15 r79342)
* macOS High Sierra 10.13.6, R 4.0.3 

## R CMD check results
There were no ERRORs or WARNINGs or NOTEs. 

## Rhub check results

  `Minor bug fixes`
  `Package has help file(s) containing install/render-stage \Sexpr{} expressions but no prebuilt PDF manual.`
  
  This is caused by the mathjaxr package which shows Latex equations in IFAA() and MZILN() documentation.
  
  `Maintainer: 'Zhigang Li <zhigang.li@ufl.edu>'`
  
## Downstream dependencies
There are currently no downstream dependencies for this package
