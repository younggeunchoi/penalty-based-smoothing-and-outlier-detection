
# Fused lasso regression with general structure

# y : Binary response vector
# X1 : Design matrix corresponding to non-regularized coefficients
# X2 : Design matrix corresponding to regularized coefficients with fused lasso penalty
# E : Set of edges that we want to give a penalty, m x 2 matrix with pairs (i,j)
# lambda1 : penalty on l_1 norm of coefficients of X2
# lambda2 : penalty on l_1 norm of differences of pairs of coefficients defined in E


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
  
  param = Xty/xsq
  param_old = param
  
  XTX = t(X)%*%X
  #if (verbose) print(XTX[1:3, 1:3])  

for(iter in 1:max_iter)
{

    # corresponding to 'A' in the paper
	  diag(A) = lambda1/(abs(param_old)+tol);

    # corresponding to 'B' in the paper
		for(j in 1:m)
		{
		  # if (verbose) cat(sprintf("j: %d\n", j))
			idx1 =E[j,1]
			idx2 = E[j,2]
      wt = E[j,3] #### added by YG
			A[idx1,idx2] = -1 * (lambda2) * wt / (abs(param_old[idx2] - param_old[idx1])+tol);
			#cat(sprintf("idx1: %d\n", idx1))
			#cat(sprintf("idx2: %d\n", idx2))
			#cat(sprintf("lambda2: %.3f\n", lambda2))
			# cat(sprintf("wt: %.3f\n", wt))	
			# cat(sprintf("param1: %.8f\n", param_old[idx1]))
			# cat(sprintf("param2: %.8f\n", param_old[idx2]))
			# cat(sprintf("abs.diff: %.8f\n", abs(param_old[idx2] - param_old[idx1])))	
				
			# if (verbose) print(-1 * (lambda2) * wt / (abs(param_old[idx2] - param_old[idx1])+tol))
        ## wt added by YG
			A[idx2,idx1] = A[idx1,idx2];
			A[idx1,idx1] = A[idx1,idx1] + (lambda2) * wt / (abs(param_old[idx2] - param_old[idx1])+tol);
      # if (verbose) print(lambda2) * wt / (abs(param_old[idx2] - param_old[idx1])+tol)  
			## wt added by yG
			A[idx2,idx2] = A[idx2,idx2] + (lambda2) * wt / (abs(param_old[idx2] - param_old[idx1])+tol);
        ## wt added by yG      
		}
	  # if (verbose) print(A[1:3, 1:3])
    
    XAmat = XTX + A ;
    # if (verbose) print(XAmat[1:3, 1:3])
    if (!pcg) {
      param = solve(XAmat, Xty)
    }
    if (pcg) {
      param = pcg(XAmat, Xty, M=diag(diag(XAmat)), 
        maxiter = 1e+05, tol = 1e-06)
    }

		max_diff = max(abs(param-param_old))

		if(max_diff <= tol)
			break;

    param_old = param

}# End Iteration


res_param = param
res_param[abs(res_param)<=tol] = 0
for(j in 1:m)
{
  	idx1 =E[j,1]
    idx2 = E[j,2]
    if(abs(res_param[idx1]-res_param[idx2])<tol)
       res_param[idx1] = res_param[idx2] = (res_param[idx1]+res_param[idx2])/2
}

return(list(param=res_param,t_iter=iter,c_tol = max_diff))
}





  
  