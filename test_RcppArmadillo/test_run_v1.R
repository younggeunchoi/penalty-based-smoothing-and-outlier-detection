library(Rcpp)
library(RcppArmadillo)

setwd("/Users/cyg/Dropbox/Codes/PHINEX/")
sourceCpp("./test_RcppArmadillo/test_source.cpp")

set.seed(1)
myx = matrix(rnorm(3*3),3 )
set.seed(2)
myy = rnorm(3)

## those two should give the same answer
a3(myx, myy) 
solve(myx, myy)

