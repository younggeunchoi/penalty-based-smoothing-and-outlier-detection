##
## @file main.R
## @brief a vignette to run the proposed method 
##	and other competitive methods
## @author Young-Geun Choi
## @note the function for the proposed method
##  is saved in solvers.R


###############################################################
################ SETTING TRUE PARAMETERS ######################
###############################################################
require(Rcpp)
require(RcppArmadillo)
setwd("/Users/cyg/Dropbox/Codes/PHINEX_pub/")

Knum = 240 # number of regions
indivnum = 50 # number of individuals per region
gam = 2 # magnitude of gamma_i's
outnum = floor(Knum * 0.10) # number of outliers

repl = 1 # replication index

###############################################################
###################### DATA GENERATION ########################
###############################################################

logit = function(t) {log(t/(1-t))}
invlogit = function(t) {exp(t)/(1+exp(t))}

## block-level informations
n.regions = Knum ;
ID.region = 1:n.regions ;
n.outliers = outnum ;
# n.outliers = 0 ;
set.seed(11*repl) ;
coord.S = runif(n=n.regions, min=5, max=95) # one-dimensional coordinate, S_i
coord.S = coord.S[order(coord.S)] ;
prop.obese.true = rep(0.4, n.regions)
prop.obese.true[coord.S > 35] = 0.5 ;
prop.obese.true[coord.S > 65] = 0.6 ;
beta = logit(prop.obese.true) ;


# set outliers (random select model)
outlier.neg = rep(FALSE, n.regions) ;
outlier.pos = rep(FALSE, n.regions) ;
gamma = rep(0, n.regions) ;
if (n.outliers > 0) {
	set.seed(22*repl) ;
	outlier.ind = sample(n.regions, size=n.outliers) ; 
	set.seed(23*repl) ;
	temp = rbinom(n=n.outliers, size=1, prob=0.5) ;
	outlier.ind.pos = outlier.ind[as.logical(temp)] ;			
	outlier.pos[outlier.ind.pos] = TRUE ;
	outlier.ind.neg = outlier.ind[!as.logical(temp)] ;		
	outlier.neg[outlier.ind.neg] = TRUE ;
}
outlier = outlier.pos | outlier.neg ;
idx_true = which(outlier) ;

gamma[outlier.pos] = gam ;
gamma[outlier.neg] = -gam ;

set.seed(33*repl) ;
# n : number of subjects enrolled
n = rep(indivnum, n.regions)
n.cum = cumsum(n) ;
## DF.region: region-level table
DF.region = data.frame(ID.region, coord.S, prop.obese.true, beta, gamma, outlier.pos, outlier.neg, 
  n)

## indiv-level
# number of people in each region. vary from 30 to 100
set.seed(44*repl)

# construct covariate information
DF.indiv = NULL ;
ID.Z2.1 = c() ;
Z2.regionLev = c() ;
for (i in 1:n.regions) {
  set.seed(i*12*repl)
  cov.indiv.temp1 = rbinom(n = n[i], size = 1, prob = 0.5) ;
  set.seed(i*22*repl)
  cov.indiv.temp2 = rbinom(n = 1, size = 1, prob = 0.5) ;
  if (cov.indiv.temp2 == 1) ID.Z2.1 = c(ID.Z2.1, cov.indiv.temp2) ; 
  Z2.regionLev = c(Z2.regionLev, cov.indiv.temp2) ;
  DF.indiv.temp = data.frame(
    ID.indiv = ifelse(length(n.cum[i-1])==0, 
      1, n.cum[i-1] + 1) : n.cum[i],
    ID.region = rep(ID.region[i], n[i]),
    Z1 = cov.indiv.temp1, Z2 = cov.indiv.temp2
    )
  DF.indiv = rbind(DF.indiv, DF.indiv.temp)
}

## DF: individual-level table
DF = merge(DF.region, DF.indiv, by="ID.region")
alpha.true = c(-0.2, 0.2) ;
# DF$weight = rep(1, nrow(DF)) ;
DF$weight = runif(n=nrow(DF), min=1, max=1) ;
DF$prop.obese.true = invlogit(DF$beta + DF$gamma + 
    as.numeric(as.matrix(DF[ ,c("Z1", "Z2")]) %*% alpha.true))
set.seed(55*repl) ;
DF$obese.ind = rbinom(n=nrow(DF), size=1, prob=DF$prop.obese.true) ;
# prevent the case of n.obese == 0 or n.obese = N: resample
for (i in 1:n.regions) {
  while(sum(DF$obese.ind[DF$ID.region == i]) == 0 | 
  	sum(DF$obese.ind[DF$ID.region == i]) == length(DF$obese.ind[DF$ID.region == i])) {
  	n_i = length(DF$obese.ind[DF$ID.region == i]) ;
  	p_i = unique(DF$prop.obese.true[DF$ID.region == i]) ;
  	DF$obese.ind[DF$ID.region == i] = rbinom(n=n_i, size=1, prob=p_i) ;
  }
}
DF.region$n = tapply(X=DF$weight, INDEX=DF$ID.region, FUN=sum)
DF.region$n.obese = tapply(X=DF$obese.ind * DF$weight, INDEX=DF$ID.region, FUN=sum)

DF.region$prop.obese = DF.region$n.obese / DF.region$n ;
DF.region$prop.obese.true.marginal = 
	invlogit( Z2.regionLev * alpha.true[2] + 
	DF.region$beta + DF.region$gamma ) / 2 +
	invlogit( alpha.true[1] + Z2.regionLev * alpha.true[2] + 
	DF.region$beta + DF.region$gamma ) / 2


# load functions
source("./src/FusedLasso_MM.R")
source("./src/solvers_debug.R") ;

## construct rho_ij
## Emat: a thin-table representation of [rho_{ij}] where
##  the 1st column for i, the 2nd column for j, the 3rd column for rho_{ij}
## Dmat: nrow(Emat) x K numeric matrix, another representation of Dmat
##  for each row in Dmat, i-th column has value rho_{ij},
##  j-th column has value -rho_{ij}, other columns has zero.
sim_dist = abs(outer(DF.region$coord.S, DF.region$coord.S, `-`))

#### Construct Dmat and Emat
DEftn = function(tvec) {
  # Construct distance between blockgroups
  sim_dist = abs(outer(tvec, tvec, `-`))
  ## try local quadratic approximation
  # construct graph E (for generalized lasso)
  sim_dist_NN = t(apply(sim_dist, 2, function(t) { t[rank(t) > 4] = 0 ; return(t) })) ;
  Emat = NULL ;
  for (i in 1:length(tvec)) { for (j in 1:length(tvec)) {
    if (sim_dist_NN[i,j] > 0) {
      if (i < j) {lowind=i ; hiind=j ;} else {lowind=j ; hiind=i ;}
      ref = which(Emat[ ,1] == lowind & Emat[ ,2] == hiind)
      if (length(ref) == 1) Emat[ref, 3] = Emat[ref, 3] + 1/(sim_dist_NN[i,j]) ;
      if (length(ref) == 0) Emat = rbind(Emat, t(c(lowind, hiind, 1/(sim_dist_NN[i,j])))) ;
    }
  }}
  if( sum(Emat[,1]>Emat[,2])>0) {
       Emat[Emat[,1]>Emat[,2],1:2] <- Emat[Emat[,1]>Emat[,2],2:1]
  }
  Emat = round(Emat, 4) ;
  Emat = unique(Emat) ;
  Emat = Emat[Emat[ ,3] > 0, ] ;
  if (is.null(dim(Emat))) Emat = t(Emat) ;
  Emat[ ,3] = Emat[ ,3] / max(Emat[ ,3]) ;
  
  # equivalent penalty matrix (D %*% beta for generalized lasso)
  Dmat = matrix(0, nrow=nrow(Emat), ncol=length(tvec)) ;
  for (i in 1:nrow(Dmat)) {
    Dmat[i, Emat[i,1]] = Emat[i,3] ;
    Dmat[i, Emat[i,2]] = -Emat[i,3] ;
  }
  return(list(Dmat=Dmat, Emat=Emat)) ;
}
temp = DEftn(DF.region$coord.S) ;
Dmat = temp$Dmat ; Emat = temp$Emat ;



###############################################################
######################## RUN METHODS ##########################
###############################################################

ymat = cbind(success = DF$obese.ind, failure = 1-DF$obese.ind)
DF$ID.region = as.factor(DF$ID.region) ;

source("./src/solvers_debug.R") ;

# single fit: no active shooting + legacy MM
source("./src/FusedLasso_MM_bak.R")
time1 = Sys.time()
res_proposed1 = try(ftn.SingleFit(Zmat = as.matrix(DF[ ,c("Z1", "Z2")]),
  ID.region = DF$ID.region,
  weight = DF$weight,
  yvec = DF$obese.ind,
  tuning_fusion = 0.6,
  tuning_sparse = 0.4,
  penalty = "hard",
  MAXIT = 100,
  verbose = TRUE,
  active_shooting=FALSE,
  Dmat = Dmat,
  Emat = Emat
  )) 
print(Sys.time() - time1)
# Time difference of 9.379224 secs
print(res_proposed1$gamma[res_proposed1$gamma != 0])


# single fit: no active shooting + MMCpp
sourceCpp("./src/FusedLasso_MM_inner.cpp")
source("./src/FusedLasso_MM.R")
time1 = Sys.time()
res_proposed2 = try(ftn.SingleFit(Zmat = as.matrix(DF[ ,c("Z1", "Z2")]),
  ID.region = DF$ID.region,
  weight = DF$weight,
  yvec = DF$obese.ind,
  tuning_fusion = 0.6,
  tuning_sparse = 0.4,
  penalty = "hard",
  MAXIT = 100,
  verbose = TRUE,
  active_shooting = FALSE,
  Dmat = Dmat,
  Emat = Emat
  )) 
print(Sys.time() - time1)
# Time difference of 4.146628 secs
print(res_proposed2$gamma[res_proposed2$gamma != 0])

# single fit: active shooting + MMCpp
sourceCpp("./src/FusedLasso_MM_inner.cpp")
source("./src/FusedLasso_MM.R")
time1 = Sys.time()
res_proposed3 = try(ftn.SingleFit(Zmat = as.matrix(DF[ ,c("Z1", "Z2")]),
  ID.region = DF$ID.region,
  weight = DF$weight,
  yvec = DF$obese.ind,
  tuning_fusion = 0.6,
  tuning_sparse = 0.4,
  penalty = "hard",
  MAXIT = 100,
  verbose = TRUE,
  active_shooting = TRUE,
  Dmat = Dmat,
  Emat = Emat
  )) 
print(Sys.time() - time1)
# Time difference of 0.8885021 secs
print(res_proposed3$gamma[res_proposed3$gamma != 0])


# no active shooting
# > res_proposed1$gamma[res_proposed1$gamma != 0]
#  [1] -1.36057  1.42774  2.63598  2.04632 -1.48918  2.48765 -1.07385 -1.69331
#  [9] -2.16415 -1.86790  2.29137  1.66416  2.18245 -2.62550  2.13485 -1.68601
# [17]  1.33014  2.34198  3.13341 -1.76278

# active shooting
#res_proposed1$gamma[res_proposed1$gamma != 0]
#  [1] -1.35043  1.42448  2.62438  2.05938 -1.48347  2.49333 -1.68042 -2.17074
#  [9] -1.86159  2.28043  1.64974  2.18730 -2.61791  2.12310 -1.67323  1.32368
# [17]  2.34869  3.15256 -1.76298


# > res_proposed1$alpha
#   ZmatZ1   ZmatZ2 
# -0.05002  0.41701 
# > res_proposed1$beta
#  [1] -0.40231 -0.40231 -0.40231 -0.40231 -0.40231 -0.40231 -0.40231 -0.46376
#  [9] -0.46376 -0.46376 -0.46376 -0.46376 -0.46376 -0.46376 -0.46376  0.31337
# [17]  0.31338  0.31338  0.31338  0.31338
# > res_proposed1$gamma
#  [1]  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000  2.00011  0.00000
#  [9]  0.00000 -1.95384  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000
# [17]  0.00000  0.00000  0.00000  0.00000



# obj_gl (beta-step)

n_temp = 120
p_temp = 240
set.seed(2)
X_temp = matrix(rnorm(n_temp * p_temp), nrow=n_temp)
set.seed(3)
y_temp = rnorm(n_temp)
w_temp = rep(1, n_temp)

source("./src/FusedLasso_MM_bak.R")
time1 = Sys.time()
obj_gl1 = fused_gen2(y = y_temp,
        X = X_temp, w = w_temp,
        E = Emat, lambda1=0, lambda2 = 0.4, pcg=FALSE, verbose=TRUE) ;
print(Sys.time() - time1)
## Time difference of 5.066974 secs
head(obj_gl1$param)

sourceCpp("./src/FusedLasso_MM_inner.cpp")
source("./src/FusedLasso_MM.R")
time1 = Sys.time()
obj_gl2 = fused_gen2(y = y_temp,
        X = X_temp, w = w_temp,
        E = Emat, lambda1=0, lambda2 = 0.4, pcg=FALSE, verbose=TRUE) ;
print(Sys.time() - time1)
head(obj_gl2$param)

