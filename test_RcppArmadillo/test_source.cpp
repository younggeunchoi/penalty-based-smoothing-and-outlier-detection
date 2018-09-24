//[[Rcpp::depends(RcppArmadillo)]]
#include "RcppArmadillo.h"
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector a1(NumericVector x) {
	return(x) ;
}

// [[Rcpp::export]]
arma::mat a2 (arma::mat x) {
	return(x) ;
}

// [[Rcpp::export]]
arma::colvec a3 (arma::mat x, arma::colvec y) {
	arma::colvec res = arma::solve(x, y) ;
	return(res) ;
}

// test A-step
// [[Rcpp::export]]
arma::mat a4 (arma::colvec param_old, 
	double tol, double lambda1, int p) {

	arma::mat A(p,p) ;
	A.zeros() ;			
	int j = 0 ;

	for (j = 0 ; j < p ; j++) {
		A(j,j) = lambda1 / (fabs(param_old(j)) + tol) ;
	}

	return(A) ;
}

