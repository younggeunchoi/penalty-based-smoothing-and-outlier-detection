// [[Rcpp::depends(RcppArmadillo)]]
#include "RcppArmadillo.h"
using namespace Rcpp;


// inputs
// 	arma::mat XTX, E
// 	arma::colvec param0, Xty
// 	double lambda1, lambda2, tol
// 	int max_iter

// outputs
// 	arma::colvec param
// 	int iter
// 	double max_diff

// [[Rcpp::export]]
List FusedLasso_MM_inner (
	arma::mat E, arma::mat XTX, 
	arma::colvec Xty, arma::colvec param0, 
	double lambda1, double lambda2, double tol,
	int max_iter
) {

	// initialize iter, j, idx1, idx2
	int iter = 0, j = 0, idx1 = 0, idx2 = 0 ;
	// initialize wt, max_diff
	double wt = 1.1, max_diff = 1.1 ;
	// initialize p as the col dimension of A
	int p = XTX.n_rows ;
	int m = E.n_rows ;

	// initialize param_old and param
	arma::colvec param_old = param0 ;
	arma::colvec param = param_old ;

	// initialize arma::mat A, XAmat
	arma::mat A(p,p) ;
	arma::mat XAmat(p,p) ;
	A.zeros() ;
	XAmat.zeros() ;			


	for (iter = 0 ; iter < max_iter ; iter++) {

		// # corresponding to 'A' in the paper
		// diag(A) = lambda1/(abs(param_old)+tol);
		for (j = 0 ; j < p ; j++) {
			A(j,j) = lambda1 / (fabs(param_old(j)) + tol) ;
		}

		// # corresponding to 'B' in the paper
		// for(j in 1:m)
		// {
		// 	idx1 =E[j,1]
		// 	idx2 = E[j,2]
		//     wt = E[j,3] 
		// 	A[idx1,idx2] = -1 * (lambda2) * wt / (abs(param_old[idx2] - param_old[idx1])+tol);				
		// 	A[idx2,idx1] = A[idx1,idx2];
		// 	A[idx1,idx1] = A[idx1,idx1] + (lambda2) * wt / (abs(param_old[idx2] - param_old[idx1])+tol);
		// 	A[idx2,idx2] = A[idx2,idx2] + (lambda2) * wt / (abs(param_old[idx2] - param_old[idx1])+tol);
		// }
		for (j = 0 ; j < m ; j++) {
			idx1 = E(j,0) - 1 ; // because E stores index starts at 1;
			idx2 = E(j,1) - 1 ; // while Cpp loop uses zero-starting index
			wt = E(j,2) ;
			
			A(idx1,idx2) = -1 * lambda2 * wt / (fabs(param_old(idx2) - param_old(idx1)) + tol) ;
			A(idx2,idx1) = A(idx1,idx2) ;
			A(idx1,idx1) = A(idx1,idx1) + lambda2 * wt / (fabs(param_old(idx2) - param_old(idx1)) + tol) ;
			A(idx2,idx2) = A(idx2,idx2) + lambda2 * wt / (fabs(param_old(idx2) - param_old(idx1)) + tol) ;
		}

		// XAmat = XTX + A ;
		XAmat = XTX + A ; 

		// param = solve(XAmat, Xty)   
		param = arma::solve(XAmat, Xty) ;

		// convergence criteria
		max_diff = arma::abs(param - param_old).max() ;
		if (max_diff <= tol) { break ; }

		// update
		param_old = param ;

	} // end iter

	List res ;
	res["param"] = param ;
	res["iter"] = iter ;
	res["max_diff"] = max_diff ;
	return(res) ;
}