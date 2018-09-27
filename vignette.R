##
## @file vignette.R
## @brief a vignette to run the proposed method 
##	and other competitive methods
## @author Young-Geun Choi
## @note the function for the proposed method
##  is saved in solvers.R


###############################################################
################ SETTING TRUE PARAMETERS ######################
###############################################################

# K: scalar, # of regions
K = 20 
# set region ID
ID.region = 1:K ;
# set outlier region (3th)
is.outlier = c(F, F, T, F, F)
i.outlier = which(is.outlier)
# K_O: scalar, # of outliers
K_O = length(i.outlier)

# nvec: vector in which the i-th element contains n_i (# of individuals in the i-th region) (i=1,...,K) 
nvec = rep(50, K) 
# gam: scalar, magnitude of nonzero gamma_i's (how much the outlier regions have departed rate from the main trend)
gam = 2 # 

# replication index, used to fix seed numbers

###############################################################
###################### DATA GENERATION ########################
###############################################################

# define logit and inverse logit function for data generation
logit = function(t) {log(t/(1-t))}
invlogit = function(t) {exp(t)/(1+exp(t))}

############## specify block-level information ################


beta = logit(rep(0.4, K)) ;
gamma = rep(0, K) ;
gamma[i.outlier] = gam ;

# Graph weight construction
#  originally constructed from the (inversed) distances between the regions' locations
Emat = data.frame(region1 = 1:(K-1), region2 = 2:K, weight = rep(1,K-1))
Dmat = cbind(diag(rep(1,K-1)),0) - cbind(0, diag(rep(1,K-1))) 
# Set Dmat
# Set Emat

### DF.region: region-level table aggreating information above
DF.region = data.frame(ID.region, beta, gamma, 
  is.outlier, n=nvec)

############ specify individual-level information #############
# number of people in each region. vary from 30 to 100
set.seed(44)

# n.cum: for internal use
n.cum = cumsum(nvec) ;

# construct covariate information
# individual-level covariate
set.seed(12)
Z1 = rbinom(n = sum(DF.region$n), size = 1, prob = 0.5) ;
# region-level covariate
set.seed(22)
Z2 = rep(rbinom(n = K, size = 1, prob = 0.5), DF.region$n)

ID.region.indiv = rep(DF.region$ID.region, DF.region$n)
DF.indiv = data.frame(ID.region = ID.region.indiv, Z1 = Z1, Z2 = Z2)



## DF: individual-level table
DF = merge(DF.region, DF.indiv, by="ID.region")

alpha = c(-0.2, 0.2) ;

# DF$weight = rep(1, nrow(DF)) ;


DF$prop.true = invlogit(DF$beta + DF$gamma + 
    as.numeric(as.matrix(DF[ ,c("Z1", "Z2")]) %*% alpha))

set.seed(55) ;
DF$y = rbinom(n=nrow(DF), size=1, prob=DF$prop.true) ;

# prevent the case of n.obese == 0 or n.obese = N: resample
for (i in 1:K) {
  while(sum(DF$y[DF$ID.region == i]) == 0 | 
  	sum(DF$y[DF$ID.region == i]) == length(DF$y[DF$ID.region == i])) {
  	n_i = length(DF$y[DF$ID.region == i]) ;
  	p_i = unique(DF$prop.true[DF$ID.region == i]) ;
  	DF$y[DF$ID.region == i] = rbinom(n=n_i, size=1, prob=p_i) ;
  }
}





###############################################################
######################## RUN METHODS ##########################
###############################################################

# load necessary packages
require(Rcpp)
require(RcppArmadillo)
# load functions
setwd("/Users/cyg/Dropbox/Codes/PHINEX_pub/")
source("./src/FusedLasso_MM.R")
source("./src/solvers_debug.R") ;

# Proposed method
res.proposed = try(ftn.SingleFit(Zmat = as.matrix(DF[ ,c("Z1", "Z2")]),
  yvec = DF$y,
  ID.region = DF$ID.region,
  tuning_fusion = 0.6,
  tuning_sparse = 0.4,
  penalty = "hard",
  MAXIT = 100,
  verbose = TRUE,
  Dmat = Dmat,
  Emat = Emat
  )) 

print(res.proposed$alpha)
#   ZmatZ1   ZmatZ2 
# -0.23150  0.17013 
print(res.proposed$beta)
#  [1] -0.42898 -0.42898 -0.42899 -0.42899 -0.42899 -0.42899 -0.42899 -0.42899
#  [9] -0.42899 -0.42899 -0.42899 -0.42899 -0.42899 -0.42899 -0.42899 -0.42899
# [17] -0.42899 -0.42899 -0.42899 -0.42899
print(res.proposed$gamma)
#  [1] 0.00000 0.00000 1.68069 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000
# [10] 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000
# [19] 0.00000 0.00000
    