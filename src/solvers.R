## @file solvers.R
## @brief functions for running the proposed method
## @author Young-Geun Choi


########### AUXILIARY FUNCTIONS ###########
logit = function(t) {log(t/(1-t))}
invlogit = function(t) {exp(t)/(1+exp(t))}
## define a soft and hard penalty function
pen_soft = function(t, lam, w = rep(1, length(t))) {
  if (lam == Inf) return(0) ;
  return(w * lam * abs(t)) ;
}
pen_hard = function(t, lam, w = rep(1, length(t))) {
  if (lam == Inf) return(0) ;
  return(w * (t < lam & t > -lam) * (lam*abs(t) - t^2/2) +
      w * (t >= lam | t <= -lam) * (lam^2/2)) ;
}
## define a soft and thresholding function
thr_soft = function (tvec, thres) {
  temp = abs(tvec) - thres ;
  temp[temp < 0] = 0 ;
  return(sign(tvec) * temp) ;
}
thr_hard = function (tvec, thres) {
  temp = tvec ;
  temp[abs(temp) < thres] = 0 ;
  return(temp) ;
}


## @fn ftn.SingleFit = function(Zmat = as.matrix(DF[ ,c("Z1", "Z2")]),
  # ID.region = DF$ID.region,
  # weight = DF$weight,
  # yvec = DF$obese.ind,
  # tuning_fusion = 0.6,
  # tuning_sparse = 0.4,
  # penalty = "hard",
  # MAXIT = 100,
  # verbose = TRUE,
  # Dmat = Dmat,
  # Emat = Emat,
  # alpha0 = NULL, beta0 = NULL, gam0 = NULL,
  # offset.alpha = NULL, offset.beta = NULL, offset.gam = NULL
##
## @brief Runs a proposed procedure with a fixed tuning parameter
## @param Zmat N x p design matrix of indiv/region level covariates, [Z_{ij}, X_i]
## @param ID.region N x 1 vector of region membership
### should take values between 1 and K(the number of regions)
## NOTE :: MAKE SURE THAT (1) the rows of
##  Zmat, ID.region, weight and yvec indicate the same invidual
##  and (2) ID.region is in increasing order 
## @param weight, N x 1 vector of individual weight, w_{ij}
## @param yvec, N x 1 vector of response, y_{ij}
## @param tuning_fusion, nonnegative scalar, lambda_1
## @param tuning_sparse, nonnegative scalar,  lambda_2
## @param penalty, character, "hard" for hard penalty, "soft" for l1 penalty
## @param MAXIT, integer, maximum number of outer iteration
## @param verbose, logical, whether to show internal calculation
## @param Emat, a thin-table representation of [rho_{ij}] where
##  the 1st column for i, the 2nd column for j, the 3rd column for rho_{ij}
## @param Dmat, nrow(Emat) x K numeric matrix, another representation of Dmat
##  for each row in Dmat, i-th column has value rho_{ij},
##  j-th column has value -rho_{ij}, other columns has zero.
## @param alpha0, beta0 and gam0, initial values for alpha, beta and gam
## @param offset.alpha, offset.beta and offset.gam, offset values (fixed if given)
##
## @return a list of results
## @return alpha, beta and gam, vectors of fitted coefficients
### (note that beta and gam is given in the order of ID.region)
## @return tuning_fusion, lambda_1 used
## @return tuning_sparse, lambad_2 used
## @return outlier, true/false indicator of the outliers
## @return df, degrees of freedom of the fit
## @return BIC, BIC value of the fit
## @return fvalueseq, for internal use
## @return tolseq, for internal use
## @return rate.est.baseline, expit(beta) 
## @return rate.est.nogam, expit(predictors + beta)
## @return rate.est.full, expit(predictors + beta + gam)
ftn.SingleFit = function(Zmat, yvec, ID.region,
  weight = rep(1,length(yvec)), 
  tuning_fusion = 0.6,
  tuning_sparse = 0.4,
  penalty = "hard",
  MAXIT = 100,
  verbose = TRUE,
  active_shooting = FALSE,
  Dmat = Dmat,
  Emat = Emat,
  alpha0 = NULL, beta0 = NULL, gam0 = NULL,
  offset.alpha = NULL, offset.beta = NULL, offset.gam = NULL) {



  # initialize parameters
  n.obese.region.w = tapply(X=yvec*weight, INDEX=ID.region, FUN=sum) ;
  n.region.w = tapply(X=weight, INDEX=ID.region, FUN=sum) ;
  prop.obese.region.w = n.obese.region.w/n.region.w ;
  n.regions = length(unique(ID.region)) ;
  if (is.null(alpha0)) alpha0 = rep(0, ncol(Zmat)) ;
  if (!is.null(offset.alpha)) alpha0 = offset.alpha ;
  if (is.null(beta0)) beta0 = logit(prop.obese.region.w) ;
  if (!is.null(offset.beta)) beta0 = offset.beta ;
  if (is.null(gam0)) {
    set.seed(1) ;
    gam0 = rnorm(n.regions) / 5 ;
  }
  if (!is.null(offset.gam)) gam0 = offset.gam ;
  alpha.old = alpha0 ;
  beta.old = beta0 ;
  gam.old = gam0 ;

  # renormalize weight: make sum to 1 (globally)
  weight.init = weight ;
  weight = weight / sum(weight) ;
  weight.region = tapply(X=weight, INDEX=ID.region, FUN=sum) ;

  # normalize penalties: sum of the weights 
  #   for the penalty should be ????
  ## (except tuning_sparse and nu)
  Dmat = Dmat / (sum(abs(Dmat)) / 2) ;
  Emat[ ,3] = Emat[ ,3] / sum(abs(Emat[ ,3])) ;

  # define loss (negative loglik) and objective function
  if (penalty == "hard") {ftn.pen = pen_hard ; ftn.thr = thr_hard ; }
  if (penalty == "soft") {ftn.pen = pen_soft ; ftn.thr = thr_soft ; }
  ftn.loss = function(alpha, beta, gam) { 
    temp.beta.expand = rep(beta, as.numeric(table(ID.region))) ;
    temp.gam.expand = rep(gam, as.numeric(table(ID.region))) ;  
    temp.lin = as.numeric(Zmat %*% alpha) + temp.beta.expand + temp.gam.expand ;
    # temp.lin = as.numeric(cbind(1,Zmat) %*% alpha) + temp.beta.expand + temp.gam.expand ;
    return( sum( 
      weight * ( log(1 + exp(temp.lin)) - yvec * temp.lin )
    ) ) ;
  }
  ftn.obj = function(alpha, beta, gam) {
    return(ftn.loss(alpha, beta, gam) + 
        tuning_fusion * sum(abs(Dmat %*% beta)) + 
        sum(ftn.pen(t=gam, lam=tuning_sparse, w=weight.region))) ;
  }
  # initial loglikelihood

  fvalue.gam.old = ftn.obj(alpha.old, beta.old, gam.old) ;
  stack.fvalue.alpha = c() ;
  stack.fvalue.beta = c() ;
  stack.fvalue.gam = c() ;
  stack.tol = c() ;

  for (iter in 1:MAXIT) {
    
    # start 'for' loop
    
    # alpha step 
    time1 = Sys.time() ;
    if (verbose) cat("alpha step.. ")
    temp.beta.expand = rep(beta.old, as.numeric(table(ID.region))) ;
    temp.gam.expand = rep(gam.old, as.numeric(table(ID.region))) ;
    if (is.null(offset.alpha)) {
      ymat = cbind(success = yvec, failure = 1-yvec) ;
      alpha.obj = glm(ymat ~ Zmat + offset(temp.offset) - 1, weight = weight, family=binomial(link="logit"),
        data=list(ymat, Zmat, weight, temp.offset = temp.beta.expand + temp.gam.expand)) ;
      alpha = coef(alpha.obj) ;
    } else { alpha = alpha0 } ;
    alpha.comp = as.numeric(Zmat %*% alpha) ;
    # alpha.comp = as.numeric(cbind(1, Zmat) %*% alpha) ;
    if (verbose) cat(sprintf("took %.2f secs\n", difftime(Sys.time(), time1, units="secs"))) ;
    fvalue.alpha = ftn.obj(alpha, beta.old, gam.old)
    stack.fvalue.alpha = c(stack.fvalue.alpha, fvalue.alpha) ;
    if (verbose) print(fvalue.alpha)
    
    # beta step
    time1 = Sys.time() ;
    if (verbose) cat("beta step.. ")
    if (is.null(offset.beta)) {
      temp.theta = alpha.comp + temp.gam.expand ;
      temp.exp = exp(temp.beta.expand + temp.theta)
      temp.a1 = weight * temp.exp / (1 + temp.exp)^2 ;
      temp.a = as.numeric(tapply(X=temp.a1, INDEX=ID.region, FUN=sum)) ;
      temp.b1 = weight * (temp.exp / (1 + temp.exp) - yvec) ;
      temp.b = as.numeric(beta.old - tapply(X=temp.b1, INDEX=ID.region, FUN=sum) / temp.a) ;
      obj_gl = fused_gen2(y = temp.b,
        X = diag(rep(1,length(temp.b))), w = temp.a,
        E = as.matrix(Emat), lambda1=0, lambda2 = tuning_fusion, pcg=FALSE, verbose=TRUE) ;
      beta = as.numeric(obj_gl$param) ;
      if (verbose) cat(sprintf("took %.2f secs\n", difftime(Sys.time(), time1, units="secs"))) ;
      fvalue.beta = ftn.obj(alpha, beta, gam.old) ;
      if (verbose) print(fvalue.beta) ;
      # modified quadratic approximation 
      if (fvalue.beta > fvalue.alpha) {
        ftnvalseq = rep(0, 100) ;
        hseq = seq(from=0.01, to=1, length=100)
        for (i in 1:length(ftnvalseq)) {
          h = hseq[i] ;
          beta.h = h * beta.old + (1-h) * beta ;
          ftnvalseq[i] = ftn.obj(alpha, beta.h, gam.old) ;
        }
        min.i = which.min(ftnvalseq) ;
        beta = hseq[min.i] * beta.old + (1-hseq[min.i]) * beta ;
        if (verbose) cat("Notice:: LQA modified\n") ;
        if (verbose) print(fvalue.beta) ;
        fvalue.beta = ftn.obj(alpha, beta, gam.old) ;
      }
    } else {
      beta = beta0 ;
      fvalue.beta = ftn.obj(alpha, beta, gam.old) ;
    }   
    stack.fvalue.beta = c(stack.fvalue.beta, fvalue.beta) ;
    
    
    # gam step
    time1 = Sys.time() ;
    if (verbose) cat("gam step.. ")    
    if (is.null(offset.gam)) {
      temp.beta.expand = rep(beta, as.numeric(table(ID.region))) ;
      temp.nu = alpha.comp + temp.beta.expand ;
      if (iter > 1) gam.prev = gam ;  
      gam = rep(0, n.regions) ;
      if (tuning_sparse != Inf) {
        for (i in 1:length(gam)) {
          if (iter > 1) { if (active_shooting & (gam.prev[i] == 0)) next ;}
          temp.nu.i = temp.nu[ID.region == i] ;
          weight.i = weight[ID.region == i] ;
          yvec.i = yvec[ID.region == i] ;
          temp.ftn.loss = function(t) { # receives vector t
            temp.lin = outer(temp.nu.i, t, "+") ;
            return( colSums( weight.i * (log(1 + exp(temp.lin)) - yvec.i * temp.lin) ) ) ;
          }
          temp.ftn.obj = function(t) {
              return(temp.ftn.loss(t) + ftn.pen(t, tuning_sparse, sum(weight.i))) ;
          }        
          
          if (penalty == "hard") {  
            gam1 = optim(fn=temp.ftn.loss, par=0, method="BFGS")$par ;
            gam2 = seq(from = -tuning_sparse, to = tuning_sparse, length=100) ;
            gam3 = 0 ;
            gams = c(gam1, gam2, gam3) ;
            gam.init = gams[which.min(temp.ftn.obj(gams))] ;           
            gam[i] = optim(fn=temp.ftn.obj, par=gam.init, method="L-BFGS-B",
              lower=-max(abs(gams)), upper=max(abs(gams)))$par ;
          }        
          if (penalty == "soft") {
            gam[i] = optim(fn=temp.ftn.obj, par=0, method="BFGS")$par ;
          }
        }
      }
    } else { gam = gam0 } ;

    if (verbose) cat(sprintf("took %.2f secs\n", difftime(Sys.time(), time1, units="secs"))) ;
    fvalue.gam = ftn.obj(alpha, beta, gam) ;
    stack.fvalue.gam = c(stack.fvalue.gam, fvalue.gam) ;
    if (verbose) print(fvalue.gam) ;
    
    # check convergence
    conv.criteria = abs(fvalue.gam.old - fvalue.gam) / pmax(abs(fvalue.gam.old), 0.1) ;
    if (verbose) { cat("tolerence: ") ; print(conv.criteria) ; }
    stack.tol = c(stack.tol, conv.criteria) ;
    if (conv.criteria < 10^-6) { break ; } else {
      alpha.old = alpha ; beta.old = beta ; gam.old = gam ;
      fvalue.gam.old = fvalue.gam ;
    }

  #end of iteration  
  }

  alpha = round(alpha, 5)
  beta = round(beta, 5) ;
  gam = round(gam, 5) ;

  ## calculate BIC

  #
  rate.est.baseline = invlogit(beta) ;
  alpha.comp = as.numeric(Zmat %*% alpha)  ;
  beta.comp = rep(beta, as.numeric(table(ID.region))) ;
  gam.comp = rep(gam, as.numeric(table(ID.region))) ;
  rate.est.nogam = tapply(
    X = weight * invlogit(alpha.comp + beta.comp),
    INDEX = ID.region, FUN = sum) /
    weight.region ;
  rate.est.full = tapply(
    X = weight * invlogit(alpha.comp + beta.comp + gam.comp),
    INDEX = ID.region, FUN = sum) /
    weight.region ;



  # calculate DF
  df.alpha = length(alpha) ;
  if (!is.null(offset.alpha)) df.alpha = 0 ;
  # beta : the number of consecutive nonzero blocks
  # df.beta = sum(beta[-1L] != beta[-length(beta)]) + 1 ;
  # in 2d, the definition should be different
  df.beta = length(unique(beta)) ;
  if (!is.null(offset.beta)) df.beta = 0 ;
  # gam : the number of nonzeros
  df.gam = sum(gam != 0)
  if (!is.null(offset.gam)) df.gam = 0 ;
  # (slightly modified) BIC: -2(loglik) + df*(log(N)+1)
  ### ref : She and Owen (2011)
  # recall : loglik() is infact a negative normalized loglik
  N = sum(weight.init) ;
  BIC = 2 * N * ftn.loss(alpha, beta, gam) + 
    (df.alpha + df.beta + df.gam) * (log(N) + 1) ;

  return(list(alpha=alpha, beta=beta, gam=gam,
    tuning_fusion=tuning_fusion, 
    tuning_sparse=tuning_sparse,
    outlier = (gam != 0),
    df = c(alpha=df.alpha, beta=df.beta, gam=df.gam),
    BIC = BIC,
    fvalueseq = stack.fvalue.gam,
    tolseq = stack.tol,
    iter=iter,
    rate.est.baseline = rate.est.baseline,
    rate.est.nogam = rate.est.nogam,
    rate.est.full = rate.est.full
    ))
}



## @fn ftn.bic = function(Zmat, yvec, ID.region,
  # weight = rep(1,length(yvec)),
  # max.nclust = Inf,
  # seq.tuning_fusion = c(2^seq(from=2, to=-5, length=100), 0),
  # seq.tuning_sparse = c(2^seq(from=2, to=-5, length=10)),
  # penalty = "hard",
  # MAXIT = 100,
  # verbose = TRUE,
  # Dmat = Dmat,
  # Emat = Emat,
  # do.plot= FALSE,
  # alpha0 = NULL, beta0 = NULL, gam0 = NULL,
  # offset.alpha = NULL, offset.beta = NULL, offset.gam = NULL) {
##
## @brief Runs a proposed procedure over 2d grid search of tuning parameters
## @param Zmat N x p design matrix of indiv/region level covariates, [Z_{ij}, X_i]
## @param ID.region N x 1 vector of region membership
### should take values between 1 and K(the number of regions)
## NOTE :: MAKE SURE THAT (1) the rows of
##  Zmat, ID.region, weight and yvec indicate the same invidual
##  and (2) ID.region is in increasing order 
## @param weight, N x 1 vector of individual weight, w_{ij}
## @param yvec, N x 1 vector of response, y_{ij}
## @param max.nclust, integer, for internal use
## @param seq.tuning_fusion, nonnegative numeric vector, sequence of lambda_1
## @param seq.tuning_sparse, nonnegative numeric vector, sequence of lambda_2
## @param penalty, character, "hard" for hard penalty, "soft" for l1 penalty
## @param MAXIT, integer, maximum number of outer iteration
## @param verbose, logical, whether to show internal calculation
## @param Emat, a thin-table representation of [rho_{ij}] where
##  the 1st column for i, the 2nd column for j, the 3rd column for rho_{ij}
## @param Dmat, nrow(Emat) x K numeric matrix, another representation of Dmat
##  for each row in Dmat, i-th column has value rho_{ij},
##  j-th column has value -rho_{ij}, other columns has zero.
## @param do.plot, for internal use
## @param alpha0, beta0 and gam0, initial values for alpha, beta and gam
## @param offset.alpha, offset.beta and offset.gam, offset values (fixed if given)
##
## @return a list of results
## @return res. the resulting object of the fit from the minimum BIC
## @return seq.tuning_fusion, the sequence of lambda_1 used
## @return seq.tuning_sparse, the sequence of lambad_2 used
## @return BIC.array, the array of of BIC values,
##  row header corresponding to seq.tuning_fusion,
##  column header corresponding to seq.tuning_sparse
## @return dfbeta.array, the array of of df(beta) values
## @return dfgam.array, the array of of df(gam) values
## @return res.array, the 2-dimensional list of fit
##  1st index corresponding to seq.tuning_fusion,
##  2nd index corresponding to seq.tuning_sparse
ftn.bic = function(Zmat, yvec, ID.region,
  weight = rep(1,length(yvec)),
  max.nclust = Inf,
  seq.tuning_fusion = c(2^seq(from=2, to=-5, length=100), 0),
  seq.tuning_sparse = c(2^seq(from=2, to=-5, length=10)),
  penalty = "hard",
  MAXIT = 100,
  verbose = TRUE,
  active_shooting = FALSE,
  Dmat = Dmat,
  Emat = Emat,
  do.plot= FALSE,
  alpha0 = NULL, beta0 = NULL, gam0 = NULL,
  offset.alpha = NULL, offset.beta = NULL, offset.gam = NULL) {

  if (penalty != "soft" & is.null(alpha0)) {
    cat("Penalty is not the soft-thresholding.\n")
    cat("Set soft-thresholding coefficients as initial values automatically.\n")
    res.init = ftn.SingleFit(Zmat = Zmat, ID.region = ID.region, 
      weight = weight, yvec = yvec, 
      tuning_fusion = seq.tuning_fusion[[1]], tuning_sparse = Inf,
      penalty = "soft", MAXIT = MAXIT, verbose = verbose,
      Dmat = Dmat, Emat = Emat,
      alpha0 = alpha0, beta0 = beta0, gam0 = gam0,
      offset.alpha = offset.alpha, offset.beta = offset.beta, offset.gam = offset.gam)
    alpha0 = res.init$alpha ; beta0 = res.init$beta ; gam0 = res.init$beta ;
  } 
  
  seq.tuning_sparse = c(Inf, seq.tuning_sparse)
  # nrow: lambda1 (fusion), ncol: lambda2 (sparse)
  array.BIC = array.dfbeta = array.dfgam = 
    matrix(0, nrow=length(seq.tuning_fusion), 
    ncol=length(seq.tuning_sparse)) 

  ### set initial objects
  # (tuning_sparse = infty and decease tuning_fusion)
  # if the number of clusters exceed certain level, then stop
  array.list.res = list(list()) ;
  
  # can mimic BIC step 1
  #list.res.step1.init = list() ;
  alpha.old = alpha0 ; beta.old = beta0 ; gam.old = gam0 ;
  # step 1
  cat("Setting tuning_sparse = infty (all gam zero)\n")
  for (k in 1:length(seq.tuning_fusion)) {
    tuning_fusion = seq.tuning_fusion[k]
    #print(tuning_fusion)
    # we need a pseudo-like dataset to get reasonable level of
    ## smoothness. go ahead for now.
    res = ftn.SingleFit(Zmat = Zmat, ID.region = ID.region, 
      weight = weight, yvec = yvec,
      tuning_fusion = tuning_fusion, tuning_sparse = Inf,
      penalty = penalty, MAXIT = MAXIT, verbose = verbose,
      active_shooting = active_shooting,
      Dmat = Dmat, Emat = Emat,
      alpha0 = alpha.old, beta0 = beta.old, gam0 = gam.old + rnorm(length(gam.old))/5,
      offset.alpha = offset.alpha, offset.beta = offset.beta, offset.gam = offset.gam) ;
    # used for initial value for next iteration
    alpha.old = res$alpha ; beta.old = res$beta ; gam.old = res$gam ;
    nclust = res$df["beta"]
    cat(sprintf("tuning_fusion=%.4f, tuning_sparse=Inf, BIC: %.3f, DF(beta,gam): %d, %d\n", 
      tuning_fusion, res$BIC, res$df["beta"], res$df["gam"])) ; 
    array.dfbeta[k,1] = res$df["beta"] ;
    array.dfgam[k,1] = res$df["gam"] ;    
    array.BIC[k,1] = res$BIC ;
    array.list.res[[k]] = list(res) ;
    if (do.plot) resplot(res, true.rate=TRUE) ;
    if (nclust > max.nclust) break ; 
  }

  ### For fixed tuning_fusion, decrease tuning_sparse
  # can mimic BIC step 2.
  # if array.dfbeta[k,1] == 0, then go next (break the loop)

  for (k in 1:length(seq.tuning_fusion)) {
    if ((array.dfbeta[k,1] == 0) & (is.null(offset.beta))) next ;
    tuning_fusion = seq.tuning_fusion[k]
    cat(sprintf("Setting tuning_fusion = %.4f\n", tuning_fusion))
    alpha.old = array.list.res[[k]][[1]]$alpha ;
    beta.old = array.list.res[[k]][[1]]$beta ;
    gam.old = array.list.res[[k]][[1]]$gam ;

    for (kk in 2:length(seq.tuning_sparse)) {
      tuning_sparse = seq.tuning_sparse[kk]
      res = ftn.SingleFit(Zmat = Zmat, ID.region = ID.region, 
        weight = weight, yvec = yvec, 
        tuning_fusion = tuning_fusion, tuning_sparse = tuning_sparse,
        penalty = penalty, MAXIT = MAXIT, verbose = verbose,
        active_shooting = active_shooting,
        Dmat = Dmat, Emat = Emat,
        alpha0 = alpha.old, beta0 = beta.old, gam0 = gam.old,
        offset.alpha = offset.alpha, offset.beta = offset.beta, offset.gam = offset.gam) ;
      # used for initial value for next iteration
      alpha.old = res$alpha ; beta.old = res$beta ; gam.old = res$gam ;
      nclust = res$df["beta"]      
      cat(sprintf("tuning_fusion=%.4f, tuning_sparse=%.4f, BIC: %.3f, DF(beta,gam): %d, %d\n", 
        tuning_fusion, tuning_sparse, res$BIC, res$df["beta"], res$df["gam"])) ; 
      array.dfbeta[k,kk] = res$df["beta"] ;
      array.dfgam[k,kk] = res$df["gam"] ;
      array.BIC[k,kk] = res$BIC ;
      array.list.res[[k]][[kk]] = res ;
      if (do.plot) resplot(res, true.rate=TRUE) ;
    } 
  }

  BIC.min = min(array.BIC) ;
  minind = which(array.BIC == BIC.min, arr.ind = TRUE)
  if (length(minind) > 2) minind = minind[1, ] ;
  res = array.list.res[[minind[1]]][[minind[2]]] ;

  return(c(res, list(
    seq.tuning_fusion = seq.tuning_fusion,
    seq.tuning_sparse = seq.tuning_sparse,
    BIC.array = array.BIC,
    dfbeta.array = array.dfbeta,
    dfgam.array = array.dfgam,
    res.array = array.list.res))) ;
}
