library(Rcpp)
library(RcppArmadillo)

setwd("/Users/cyg/Dropbox/Codes/PHINEX_pub/")
# sourceCpp("./test_Rcpp/test_matrix.cpp")

# a1(c(1,2,3))

# a2(as.matrix(c(1,2,3)));

# myx = matrix(rnorm(3*3),3 )
# myy = rnorm(3)
# a3(myx, myy)

# solve(myx, myy)

# a3(myx, myy)

# # arma::mat a4 (arma::colvec param_old, 
# # 	double tol, double lambda1, int p)
# a4(param_old = c(1,2,3), tol = 1, lambda1 = 1, p=3)
# a4(c(1,2,3), tol = 1, lambda1 = 1, p=3)


sourceCpp("./src/FusedLasso_MM_inner.cpp")

# List FusedLasso_MM_inner (
# 	arma::mat E, arma::mat XTX, 
# 	arma::colvec param0, arma::colvec Xty,
# 	double lambda1, double lambda2, double tol,
# 	int max_iter

## test at simple settings
p = 5 ;
E = rbind(c(1,2,1), c(2,3,1), c(3,4,1), c(4,5,1)) ;
XTX = diag(rep(1,p)) ;
set.seed(1) ;
Xty = rnorm(p) ;
param0 = rep(0,p) ;
lambda1 = 0.1 ; lambda2 = 0.2 ;
tol = 10^-7 ;
max_iter = 1000 ;

time1 = Sys.time() ;
res = FusedLasso_MM_inner (E, XTX, 
 	Xty, param0, 
 	lambda1, lambda2, tol,
 	max_iter)
print(Sys.time() - time1)

# to compare


time1 = Sys.time()
param_old = param0
A = matrix(0, p, p)
for(iter in 1:max_iter)
{

    # corresponding to 'A' in the paper
	  diag(A) = lambda1/(abs(param_old)+tol);

    # corresponding to 'B' in the paper
		for(j in 1:nrow(E))
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
    param = solve(XAmat, Xty)
    
		max_diff = max(abs(param-param_old))

		if(max_diff <= tol)
			break;

    param_old = param

}# End Iteration
print(Sys.time() - time1)
print(param)