##
## @file vignette.R
## @brief a vignette to run the proposed method 
##	and other competitive methods
## @author Young-Geun Choi
## @note the function for the proposed method
##  is saved in solvers.R


# modify the working directory
setwd("/Users/cyg/Dropbox/Codes/PHINEX_pub/")

# define logit and inverse logit function for data generation
logit = function(t) {log(t/(1-t))}
invlogit = function(t) {exp(t)/(1+exp(t))}


###############################################################
################ SETTING TRUE PARAMETERS ######################
###############################################################

# K: scalar, # of regions
K = 20 
# set region ID
ID.region = 1:K ;
# set outlier region (3th)
# K_O: scalar, # of outliers
i.outlier = 3 ;
is.outlier = rep(FALSE, K) ; is.outlier[i.outlier] = TRUE ;
K_O = length(i.outlier) ;

# nvec: vector in which the i-th element contains n_i (# of individuals in the i-th region) (i=1,...,K) 
nvec = rep(50, K) 



###############################################################
###################### DATA GENERATION ########################
###############################################################

############## specify region-level information ###############

# beta: vector of length K, specifies the true \beta_i's
beta = logit(rep(0.4, K)) ;
# gam: vector of length K, specifies the true \gamma_i's 
## (named as 'gam' to avoid the conflict with the R-builtin function gamma())
gam = rep(0, K) ;
gam[i.outlier] = 2 ;

### DF.region: region-level table aggreating information above
DF.region = data.frame(ID.region, beta, gam, 
  is.outlier, n=nvec)

############ specify individual-level information #############

# construct covariate information
# individual-level covariate
set.seed(12)
Z1 = rbinom(n = sum(DF.region$n), size = 1, prob = 0.5) ;
# region-level covariate
set.seed(22)
Z2 = rep(rbinom(n = K, size = 1, prob = 0.5), DF.region$n)
ID.region.indiv = rep(DF.region$ID.region, DF.region$n)
DF.indiv = data.frame(ID.region = ID.region.indiv, Z1 = Z1, Z2 = Z2)

####### merge regional- and indiv-level tables ################
DF = merge(DF.region, DF.indiv, by="ID.region")

# model p_ij = Pr(Y=1 | Z_ij, X_i)
alpha = c(-0.2, 0.2) ;
DF$prop.true = invlogit(DF$beta + DF$gam + 
    as.numeric(as.matrix(DF[ ,c("Z1", "Z2")]) %*% alpha))

# generate y_ij ~ Bernoulli(p_ij)
set.seed(55) ;
DF$y = rbinom(n=nrow(DF), size=1, prob=DF$prop.true) ;

# prevent the existance of i-th region such that [y_ij == 0 or == 1 for all j]
for (i in 1:K) {
  while(sum(DF$y[DF$ID.region == i]) == 0 | 
  	sum(DF$y[DF$ID.region == i]) == length(DF$y[DF$ID.region == i])) {
  	n_i = length(DF$y[DF$ID.region == i]) ;
  	p_i = unique(DF$prop.true[DF$ID.region == i]) ;
  	DF$y[DF$ID.region == i] = rbinom(n=n_i, size=1, prob=p_i) ;
  }
}

######### model spatial depndence over regions ################
#  originally constructed from the (inversed) distances between the regions' locations
Emat = data.frame(region1 = 1:(K-1), region2 = 2:K, weight = rep(1,K-1))
Dmat = cbind(diag(rep(1,K-1)),0) - cbind(0, diag(rep(1,K-1))) 
# Set Dmat
# Set Emat



###############################################################
######################## RUN METHODS ##########################
###############################################################

# load necessary packages
require(Rcpp)
require(RcppArmadillo)
# load functions
sourceCpp("./src/FusedLasso_MM_inner.cpp")
source("./src/FusedLasso_MM.R")
source("./src/solvers.R") ;

# Proposed method with a fixed tuning parameter (lambda_1 = 0.6, lambda_2 = 0.4)
res.proposed = ftn.SingleFit(Zmat = as.matrix(DF[ ,c("Z1", "Z2")]),
  yvec = DF$y,
  ID.region = DF$ID.region,
  tuning_fusion = 0.3,
  tuning_sparse = 0.3,
  penalty = "hard",
  MAXIT = 100,
  verbose = TRUE,
  Dmat = Dmat,
  Emat = Emat
  ) 

print(res.proposed$alpha)
#   ZmatZ1   ZmatZ2 
# -0.22379  0.26501 
print(res.proposed$beta)
#  [1] -0.43233 -0.43233 -0.43233 -0.43233 -0.43233 -0.43233 -0.43238 -0.43241
#  [9] -0.43241 -0.43242 -0.43242 -0.43242 -0.43242 -0.43242 -0.43242 -0.43242
# [17] -0.43242 -0.43242 -0.43242 -0.43242
print(res.proposed$gam)
#  [1]  0.00000  0.00000  1.58384  0.00000  0.00000  0.00000 -0.98093  0.00000
#  [9]  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000
# [17]  0.00000  0.00000  0.00000  0.00000
    

# Proposed method with a tuning parameter search (criteria: modified BIC as defined in the manuscript)
res.proposed.bic = ftn.bic(Zmat = as.matrix(DF[ ,c("Z1", "Z2")]),
  yvec = DF$y,
  ID.region = DF$ID.region,
  seq.tuning_fusion = c(2^seq(from=2, to=-5, length=100), 0),
  seq.tuning_sparse = c(2^seq(from=2, to=-5, length=10)),
  penalty = "hard",
  MAXIT = 100,
  verbose = FALSE,
  Dmat = Dmat,
  Emat = Emat
  ) 

print(res.proposed.bic$alpha)
#   ZmatZ1   ZmatZ2 
# -0.22714  0.17491
print(res.proposed.bic$beta)
#  [1] -0.43329 -0.43329 -0.43329 -0.43329 -0.43329 -0.43329 -0.43329 -0.43329
#  [9] -0.43329 -0.43329 -0.43329 -0.43329 -0.43329 -0.43329 -0.43329 -0.43329
# [17] -0.43329 -0.43329 -0.43329 -0.43329
print(res.proposed.bic$gam)
#  [1] 0.00000 0.00000 1.67721 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000
# [10] 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000
# [19] 0.00000 0.00000

