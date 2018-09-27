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

# load necessary packages
require(Rcpp)
require(RcppArmadillo)
setwd("/Users/cyg/Dropbox/Codes/PHINEX_pub/")

# K: scalar, # of regions
K = 240 
# nvec: vector in which the i-th element contains n_i (# of individuals in the i-th region) (i=1,...,K) 
nvec = rep(50, K) 
# gam: scalar, magnitude of nonzero gamma_i's (how much the outlier regions have departed rate from the main trend)
gam = 2 # 
# K_O: scalar, # of outliers
K_O = floor(K * 0.10) 

# replication index, used to fix seed numbers
repl = 1 

###############################################################
###################### DATA GENERATION ########################
###############################################################

# define logit and inverse logit function for data generation
logit = function(t) {log(t/(1-t))}
invlogit = function(t) {exp(t)/(1+exp(t))}

############## specify block-level information ################

# set region ID
ID.region = 1:K ;

# define S_i, geographic location to the i-th region (i=1,...,K)
# assume S_i \in [0, 100] (all the regions lies at one-dimensional line)
# Svec:= (S_1, ..., S_K)
set.seed(11*repl) ;
Svec = runif(n=K, min=5, max=95) 
# give lower S_i to lower i
Svec = Svec[order(Svec)] ;

# assign beta as follows:
#  if S_i < 35,       then beta_i = logit(0.4)
#  if 35 <= S_i < 65, then beta_i = logit(0.5)
#  if S_i >= 65,      then beta_i = logit(0.6)
beta = logit(rep(0.4, K))
beta[Svec > 35] = logit(0.5) ;
beta[Svec > 65] = logit(0.6) ;

# assign gamma as follows:
#  1) assign gamma_i =  gam to approx. [K_O/2] of i's
#  2) assign gamma_i = -gam to approx. K_O - [K_O/2] of i's
# in code below,
#  i.outlier: vector containing indices of outlier regions (gamma_i != 0)
#  i.pos:     vector containing indices of gamma_i > 0
#  i.neg:     vector containing indices of gamma_i < 0
#  is.outlier, is.pos and is.neg: logical indices of i.outlier, i.pos and i.neg
gamma = rep(0, K) ;
is.neg = rep(FALSE, K) ;
is.pos = rep(FALSE, K) ;
if (K_O > 0) {
	set.seed(22*repl) ;
  # sample a set of i's of size K_O from {1,...,K}
	i.outlier = sample(K, size=K_O) ; 
	set.seed(23*repl) ;
  # further sample a subset of size [K_O/2] that will be assigned to positive gamma
	temp = rbinom(n=K_O, size=1, prob=0.5) ;
	i.pos = i.outlier[as.logical(temp)] ;			
	is.pos[i.pos] = TRUE ;
	i.neg = i.outlier[!as.logical(temp)] ;		
	is.neg[i.neg] = TRUE ;
}
if (K_O == 0) i.outlier = c() ;
gamma[is.pos] = gam ;
gamma[is.neg] = -gam ;
is.outlier = is.pos | is.neg ;



## DF.region: region-level table
DF.region = data.frame(ID.region, Svec, beta, gamma, is.pos, is.neg, 
  n=nvec)

############ specify individual-level information #############
# number of people in each region. vary from 30 to 100
set.seed(44*repl)

# construct covariate information
DF.indiv = NULL ;
ID.Z2.1 = c() ;
Z2.regionLev = c() ;
for (i in 1:K) {
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
DF$prop.true = invlogit(DF$beta + DF$gamma + 
    as.numeric(as.matrix(DF[ ,c("Z1", "Z2")]) %*% alpha.true))
set.seed(55*repl) ;
DF$obese.ind = rbinom(n=nrow(DF), size=1, prob=DF$prop.true) ;
# prevent the case of n.obese == 0 or n.obese = N: resample
for (i in 1:K) {
  while(sum(DF$obese.ind[DF$ID.region == i]) == 0 | 
  	sum(DF$obese.ind[DF$ID.region == i]) == length(DF$obese.ind[DF$ID.region == i])) {
  	n_i = length(DF$obese.ind[DF$ID.region == i]) ;
  	p_i = unique(DF$prop.true[DF$ID.region == i]) ;
  	DF$obese.ind[DF$ID.region == i] = rbinom(n=n_i, size=1, prob=p_i) ;
  }
}
DF.region$n = tapply(X=DF$weight, INDEX=DF$ID.region, FUN=sum)
DF.region$n.obese = tapply(X=DF$obese.ind * DF$weight, INDEX=DF$ID.region, FUN=sum)

DF.region$prop = DF.region$n.obese / DF.region$n ;
DF.region$prop.true.marginal = 
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
sim_dist = abs(outer(DF.region$Svec, DF.region$Svec, `-`))

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
temp = DEftn(DF.region$Svec) ;
Dmat = temp$Dmat ; Emat = temp$Emat ;




###############################################################
######################## RUN METHODS ##########################
###############################################################

ymat = cbind(success = DF$obese.ind, failure = 1-DF$obese.ind)
DF$ID.region = as.factor(DF$ID.region) ;


# GLMM
require(lme4)
res.glmm = glmer(ymat ~ 1 + Z1 + Z2 + (1|ID.region),
  family=binomial, 
  weights=weight, 
  data=list(ymat=ymat, Z1=DF$Z1, Z2=DF$Z2, weight=DF$weight,
    ID.region = DF$ID.region))  


# Scan statistic, tweaking the algorithm from SpatialEpi
require(Rcpp)
sourceCpp("./SpatialEpi_Adj/cluster_detection_YGadj.cpp") # 
source("./SpatialEpi_Adj/pkg_adj.R") 
obj_scan = kulldorff(
  geo = cbind(DF.region$Svec, 0),
  cases = DF.region$n.obese, population = DF.region$n,
  expected.cases = DF.region$n * (sum(DF.region$n.obese)/sum(DF.region$n)), 
  pop.upper.bound = 0.5,
  n.simulations = 9999, alpha.level = 0.05, plot=FALSE)

# Proposed method
res_proposed = try(ftn.2d(Zmat = as.matrix(DF[ ,c("Z1", "Z2")]),
  ID.region = DF$ID.region,
  weight = DF$weight,
  yvec = DF$obese.ind,
  max.nclust = Inf,
  seq.tuning_fusion = c(2^seq(from=12, to=-2, length=50), 0),
  seq.tuning_sparse = c(2^seq(from=2, to=-5, length=50)),
  penalty = "hard",
  MAXIT = 100,
  verbose = FALSE,
  Dmat = Dmat,
  Emat = Emat
  )) 


    
###############################################################
######################## EVALUATION ###########################
###############################################################

# evaluate outliers
ftn_eval = function(esxvec, truevec = idx_true, all = DF.region$ID.region) {
  #(TPR = TP/P, TNR= TN/N, FDR= FP/(TP+FP), accuracy=(TP+TN)/(P+N), 
  #MCC = TP*TN - FP*FN / sqrt((TP+FP)(TP+FN)(TN+FP)(TN+FN))
  Pvec = truevec ; P = length(Pvec) ;
  Nvec = setdiff(all, truevec) ; N = length(Nvec) ;
  TP = length(intersect(Pvec, esxvec))
  FN = length(setdiff(Pvec, esxvec))
  TN = length(intersect(Nvec, setdiff(all, esxvec)))
  FP = length(setdiff(Nvec, setdiff(all, esxvec)))
  
  TPR = TP / P ; TNR = TN / N ; FDR = FP/(TP+FP) ;
  FOR = FN/(TN+FN) ;
  # FDR MAY BE ZERO
  ACC = (TP + TN) / (P + N) ;
  MCC = (TP*TN - FP*FN) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)) ;
  return(c(TPR = TPR, TNR = TNR, FDR = FDR, FOR = FOR,
    ACC = ACC, MCC = MCC)) ;
}

# evaluate coefficients
ftn_eval2 = function(est_alpha, est_beta, est_prev=NA,
  true_alpha = alpha.true, true_beta = beta, true_prev = DF.region$prop.true.marginal) {
  return(c(bias_alpha1 = as.numeric(est_alpha[1] - true_alpha[1]),
    bias_alpha2 = as.numeric(est_alpha[2] - true_alpha[2]),
    bias_beta = mean(est_beta - true_beta),
    MSE_beta = mean((est_beta - true_beta)^2),
    MSE_prev = mean((est_prev - true_prev)^2)
    ))
} 


# summarizing the results

# results from the proposed method
if (class(res_proposed) != "try-error") { 
  idx_proposed = which(res_proposed$is.outlier)
  tbl_proposed = data.frame(est="proposed", t(ftn_eval(idx_proposed)), 
    t(ftn_eval2(res_proposed$alpha, res_proposed$beta, res_proposed$rate.est.full)),
    df.beta = res_proposed$df[2]) ;
} else tbl_proposed = NULL ;

# results from the GLMM
glmm_ranef = as.numeric(unlist(ranef(res.glmm)[[1]]))
glmm_fitbeta = as.numeric(fixef(res.glmm)[1]) + glmm_ranef
glmm_sd = sqrt(as.numeric((VarCorr(res.glmm)[[1]])))
glmm_outind = which( abs(glmm_ranef / glmm_sd)  > 2.5)
glmm_est_prev = tapply(X=predict(res.glmm, type="response") * DF$weight, 
  INDEX=DF$ID.region, FUN=sum) / tapply(X=DF$weight, 
  INDEX=DF$ID.region, FUN=sum)
tbl_glmm = data.frame(est="glmm", t(ftn_eval(glmm_outind)), 
  t(ftn_eval2(fixef(res.glmm)[-1], glmm_fitbeta, glmm_est_prev)), df.beta=NA) ;

# results from the scan statistic
idx_scan = obj_scan$most.likely.cluster$location.IDs.included
tbl_scan = data.frame(est="scan", t(ftn_eval(idx_scan)), 
  bias_alpha1=NA, bias_alpha2=NA, bias_beta=NA, MSE_beta=NA, MSE_prev=NA, df.beta=NA);
    

tbl = data.frame(K_O=K_O,
  Knum = Knum, indivnum=indivnum, gam=gam, 
  rbind(tbl_scan, tbl_glmm, tbl_proposed))

print(tbl)
