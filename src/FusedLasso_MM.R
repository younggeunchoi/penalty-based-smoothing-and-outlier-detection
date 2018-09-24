
# Fused lasso regression with general structure

# y : Binary response vector
# X1 : Design matrix corresponding to non-regularized coefficients
# X2 : Design matrix corresponding to regularized coefficients with fused lasso penalty
# E : Set of edges that we want to give a penalty, m x 2 matrix with pairs (i,j)
# lambda1 : penalty on l_1 norm of coefficients of X2
# lambda2 : penalty on l_1 norm of differences of pairs of coefficients defined in E


# require to load FusedLasso_MM_inner.cpp
require(RcppArmadillo)
require(Rcpp)
sourceCpp("./src/FusedLasso_MM_inner.cpp")

fused_gen2 <- function(y,X,w = rep(1,length(y)),
  E,lambda1,lambda2,tol=1e-6,max_iter = 1000,
  pcg=FALSE, verbose=FALSE)
{
  # minimizes \sum w_i * (Y_i - X_i^T\beta)^2 + 
  # Penalty(beta;lambda1,lambda2).
  if(pcg) require(pcg) ;

  n <- length(y)
  p <- ncol(X)
  X <- sqrt(w) * X ;
  y <- sqrt(w) * y ;
   
  A = matrix(0,p,p)
	Ainv = matrix(0,p,p)

  # increasing order in each row
  if (sum(E[,1]>E[,2]) > 0)
     E[E[,1]>E[,2],1:2] <- E[E[,1]>E[,2],2:1]

  E <- unique(E)  ##

  m <- nrow(E) # number of unique edges #

  Xty = t(X) %*% y
  #print(t(X)[1:3, 1:3])
  #print(y)
  #print(head(y))
  #print(head(Xty))

  xsq = apply(X^2,2,sum)
  #print(head(xsq))
  
  param0 = Xty/xsq
  
  XTX = t(X)%*%X
  #if (verbose) print(XTX[1:3, 1:3])  


  res_inner = FusedLasso_MM_inner (E, XTX, 
    Xty, param0, lambda1, lambda2, tol, max_iter)

  res_param = res_inner$param
  res_param[abs(res_param)<=tol] = 0
  for (j in 1:m)
  {
  	idx1 =E[j,1]
    idx2 = E[j,2]
    if(abs(res_param[idx1]-res_param[idx2])<tol)
       res_param[idx1] = res_param[idx2] = (res_param[idx1]+res_param[idx2])/2
  }

return(list(param=res_param,t_iter=res_inner$iter,c_tol = res_inner$max_diff))
}





  
  