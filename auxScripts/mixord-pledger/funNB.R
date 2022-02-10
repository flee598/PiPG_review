# funNB.R 

# Functions for count data, pattern detection models.
# Negative binomial version.
# Variance modelled as a quadratic in mean.
# Basic parameters are lin. model for mu(ij) and phi(j)
# where var = mu + phi*mu^2

# xlogx function:
# ---------------

# Function for x.log(x) with zero if x=0 (or if computer puts x
# just below zero). 
# Argument may be a scalar, vector or matrix.

xlogx.fn <- function(x){
  if (is.matrix(x)) xlx <- matrix(0,nrow(x),ncol(x))
  if (is.vector(x)) xlx <- rep(0,length(x))
  xlx[x>0] <- x[x>0] * log(x[x>0])
  xlx
}


# ----------- MODEL FITTING FUNCTIONS -----------------------------

# ------------------------------------------------------------------
#               GLMs
# ------------------------------------------------------------------


# ------------ Null Model ---------------------------------------------

# GLM with log mu(ij) = lambda (scalar).
# phi(j) = coeff. in quadratic for column j.
# sigsq(ij) = mu(ij) + phi(j)*mu(ij)^2
# k = 1/phi
# Linear predictor: log(mu(ij)) = lambda

# Input vector for optim:
# lambda value  [1]
# log(k) vector [2:(p+1)]

NB.Null.ll <- function(invect){
  # Read in parameters:
  lambda.s <- invect[1]
  k.v     <- exp(invect[2:(p+1)])
  k.m     <- matrix(k.v,n,p,byrow=T)
  # Find the linear predictor and mu:
  lin.m   <- matrix(lambda.s,n,p)
  mu.m    <- exp(lin.m)
  # Calculate NB likelihood:
  # ll <- -1000000
  # if (min(c(k.m,mu.m))>0) 
  ll <-
    sum(lgamma(y.mat+k.m) - lgamma(k.m) - lgamma(y.mat+1) +
        k.m*log(k.m) + y.mat*log(mu.m) -
        (k.m + y.mat)*log(k.m + mu.m))
  ll
}

NB.Null.fn <- function(start.v=rep(0,p+1)){
  # Run the optimisation, with repeats if necessary:
  this.fit <- optim(par=start.v,
                    fn=NB.Null.ll,
                    method = "L-BFGS-B",
                    control = list(fnscale=-1))
  print(paste("Fit",1," con =",this.fit$con))
  for (i in 2:10) if (this.fit$con != 0){
    start.vect <- this.fit$par
    this.fit <- optim(par=start.vect,
                      fn=NB.Null.ll,
                      method = "L-BFGS-B",
                    control = list(fnscale=-1))
    print(paste("Fit",i," con =",this.fit$con))
  }
  # Save results:
  max.ll  <- this.fit$value
  res.dev <- -2*max.ll
  npar    <- length(this.fit$par)
  aic  <- res.dev + 2*npar
  aicc <- NA
  if (n*p > npar + 2)
    aicc <- aic + (2*(npar+1)*(npar+2))/(n*p - npar - 2)
  bic  <- res.dev + npar*log(n*p)
  icl  <- bic
  outvect    <- this.fit$par
  lambda.hat <- outvect[1]
  k.hat      <- exp(outvect[2:(p+1)])
  lin.hat    <- matrix(lambda.hat,n,p)
  mu.hat     <- exp(lin.hat)
  phi.hat    <- 1/k.hat
  sigsq.hat  <- mu.hat + mu.hat^2*matrix(phi.hat,n,p,byrow=T)
  k.hat[k.hat>1000] <- 1/0
  out1 <- round(c(n,p,max.ll,res.dev,npar,aic,aicc,bic,icl,max.ll),3)
  names(out1) <- c("n","p","Max.ll","Res.Dev.","npar","AIC","AICc",
                   "BIC","ICL","LLc")
  this.out <- list("info"=out1,
                   "lambda"=round(lambda.hat,3),
                   "phi"=round(phi.hat,3),
                   "k"=round(k.hat,3),
                   "mu"=round(mu.hat[,1],3),
                   "sigma-sq"=round(sigsq.hat[1,],3))
  this.out
}


# ------------ Model A ---------------------------------------------

# GLM with log mu(ij) = alpha(i).
# phi(j) = coeff. in quadratic for column j.
# sigsq(ij) = mu(ij) + mu(ij)^2*phi(j)

# Linear predictor: log(mu(ij)) = alpha(i)

# Input vector for optim:
# alpha vector  [1:n]
# log(k) vector [(n+1):(n+p)]

NB.A.ll <- function(invect){
  # Read in parameters:
  alpha.v <- invect[1:n]
  k.v     <- exp(invect[(n+1):(n+p)])
  k.m     <- matrix(k.v,n,p,byrow=T)
  # Find the linear predictor and mu:
  lin.m   <- matrix(alpha.v,n,p)
  mu.m    <- exp(lin.m)
  # Calculate NB likelihood:
  # ll <- -1000000
  # if (min(c(k.m,mu.m))>0) 
  ll <-
    sum(lgamma(y.mat+k.m) - lgamma(k.m) - lgamma(y.mat+1) +
        k.m*log(k.m) + y.mat*log(mu.m) -
        (k.m + y.mat)*log(k.m + mu.m))
  ll
}

NB.A.fn <- function(start.v=rep(0,n+p)){
  # Run the optimisation, with repeats if necessary:
  this.fit <- optim(par=start.v,
                    fn=NB.A.ll,
                    method = "L-BFGS-B",
                    control = list(fnscale=-1))
  print(paste("Fit",1," con =",this.fit$con))
  for (i in 2:10) if (this.fit$con != 0){
    start.vect <- this.fit$par
    this.fit <- optim(par=start.vect,
                      fn=NB.A.ll,
                      method = "L-BFGS-B",
                    control = list(fnscale=-1))
    print(paste("Fit",i," con =",this.fit$con))
  }
  # Save results:
  max.ll  <- this.fit$value
  res.dev <- -2*max.ll
  npar    <- length(this.fit$par)
  aic  <- res.dev + 2*npar
  aicc <- NA
  if (n*p > npar + 2)
    aicc <- aic + (2*(npar+1)*(npar+2))/(n*p - npar - 2)
  bic  <- res.dev + npar*log(n*p)
  icl  <- bic
  outvect   <- this.fit$par
  alpha.hat <- outvect[1:n]
  k.hat     <- exp(outvect[(n+1):(n+p)])
  lin.hat   <- matrix(alpha.hat,n,p)
  mu.hat    <- exp(lin.hat)
  phi.hat   <- 1/k.hat
  sigsq.hat <- mu.hat + mu.hat^2*matrix(phi.hat,n,p,byrow=T)
  k.hat[k.hat>1000] <- 1/0
  out1 <- round(c(n,p,max.ll,res.dev,npar,aic,aicc,bic,icl,max.ll),3)
  names(out1) <- c("n","p","Max.ll","Res.Dev.","npar","AIC","AICc",
                   "BIC","ICL","LLc")
  this.out <- list("info"=out1,
                   "alpha"=round(alpha.hat,3),
                   "phi"=round(phi.hat,3),
                   "k"=round(k.hat,3),
                   "mu"=round(mu.hat[,1],3),
                   "sigma-sq"=round(sigsq.hat[1,],3))
  this.out
}


# ------------ Model B ---------------------------------------------

# GLM with only column effects in mu.
# log mu(ij) = beta(j).
# phi(j) = coeff. in quadratic for column j.
# sigsq(ij) = mu(ij) + mu(ij)^2*phi(j)

# Linear predictor: log(mu(ij)) = alpha(i)

# Input vector for optim:
# beta vector   [1:p]
# log(k) vector [(p+1):(2*p)]

NB.B.ll <- function(invect){
  # Read in parameters:
  beta.v  <- invect[1:p]
  k.v     <- exp(invect[(p+1):(2*p)])
  k.m     <- matrix(k.v,n,p,byrow=T)
  # Find the linear predictor:
  lin.m   <- matrix(beta.v,n,p,byrow=T)
  # Find mu:
  mu.m     <- exp(lin.m)
  # Calculate NB likelihood:
  ll <- -1000000
  if (min(c(k.m,mu.m))>0) ll <-
    sum(lgamma(y.mat+k.m) - lgamma(k.m) - lgamma(y.mat+1) +
        k.m*log(k.m) + y.mat*log(mu.m) -
        (k.m + y.mat)*log(k.m + mu.m))
  ll
}

NB.B.fn <- function(start.v=rep(0,2*p)){
  # Run the optimisation, with repeats if necessary:
  this.fit <- optim(par=start.v,
                    fn=NB.B.ll,
                    method = "L-BFGS-B",
                    control = list(fnscale=-1))
  print(paste("Fit",1," con =",this.fit$con))
  for (i in 2:10) if (this.fit$con != 0){
    start.vect <- this.fit$par
    this.fit <- optim(par=start.vect,
                      fn=NB.B.ll,
                      method = "L-BFGS-B",
                      control = list(fnscale=-1))
    print(paste("Fit",i," con =",this.fit$con))
  }
  # Save results:
  max.ll  <- this.fit$value
  res.dev <- -2*max.ll
  npar    <- length(this.fit$par)
  aic  <- res.dev + 2*npar
  aicc <- NA
  if (n*p > npar + 2)
    aicc <- aic + (2*(npar+1)*(npar+2))/(n*p - npar - 2)
  bic  <- res.dev + npar*log(n*p)
  icl  <- bic
  outvect   <- this.fit$par
  beta.hat  <- outvect[1:p]
  k.hat     <- exp(outvect[(p+1):(2*p)])
  phi.hat   <- 1/k.hat
  lin.hat   <- matrix(beta.hat,n,p,byrow=T)
  mu.hat    <- exp(lin.hat)
  sigsq.hat <- mu.hat + mu.hat^2*matrix(phi.hat,n,p,byrow=T)
  k.hat[k.hat>1000] <- 1/0
  out1 <- round(c(n,p,max.ll,res.dev,npar,aic,aicc,bic,icl,max.ll),3)
  names(out1) <- c("n","p","Max.ll","Res.Dev.","npar","AIC","AICc",
                   "BIC","ICL","LLc")
  this.out <- list("info"=out1,
                   "beta"=round(beta.hat,3),
                   "phi"=round(phi.hat,3),
                   "k"=round(k.hat,3),
                   "mu"=round(mu.hat[1,],3),
                   "sigma-sq"=round(sigsq.hat[1,],3))
  this.out
}


# ------------ Model AB -----------------------------------------

# Basic null model of no association.
# log mu(ij) = alpha(i) + beta(j) with constraint sum alpha(i)=0.
# phi(j) = coeff. in quadratic for column j.
# sigsq(ij) = mu(ij) + mu(ij)^2*phi(j)
# theta(ij) = mu(ij)/sigsq(ij)
# eta(ij) = mu(ij)*theta(ij)/(1-theta(ij))

# Parameters for optim:
# alpha[-n]   [1:(n-1)]
# beta        [n:(n+p-1)]
# log(k)      [(n+p):(n+2*p-1)]

NB.AB.ll <- function(invect){
  # Read in parameters:
  alpha.v <- invect[1:(n-1)]
  alpha.v <- c(alpha.v,0-sum(alpha.v))
  beta.v  <- invect[n:(n+p-1)]
  k.v     <- exp(invect[(n+p):(n+2*p-1)])
  k.m     <- matrix(k.v,n,p,byrow=T)
  phi.m   <- 1/k.m
  # Find the linear predictor:
  lin.m   <- outer(alpha.v,beta.v,FUN="+")
  # Find mu and sigsq:
  mu.m     <- exp(lin.m)
  sigsq.m  <- mu.m + mu.m^2*phi.m
  # Calculate NB likelihood:
  ll <- -1000000
  if (min(c(k.m,mu.m))>0) ll <-
    sum(lgamma(y.mat+k.m) - lgamma(k.m) - lgamma(y.mat+1) +
        k.m*log(k.m) + y.mat*log(mu.m) -
        (k.m + y.mat)*log(k.m + mu.m))
  ll
}

NB.AB.fn <- function(start.v=rep(0,n+2*p-1)){
  # Run the optimisation, with repeats if necessary:
  this.fit <- optim(par=start.v,
                    fn=NB.AB.ll,
                    method = "L-BFGS-B",
                    control = list(fnscale=-1))
  print(paste("Fit",1," con =",this.fit$con))
  for (i in 2:10) if (this.fit$con != 0){
    start.vect <- this.fit$par
    this.fit <- optim(par=start.vect,
                      fn=NB.AB.ll,
                      method = "L-BFGS-B",
                    control = list(fnscale=-1))
    print(paste("Fit",i," con =",this.fit$con))
  }
  # Save results:
  max.ll  <- this.fit$value
  res.dev <- -2*max.ll
  npar    <- length(this.fit$par)
  aic  <- res.dev + 2*npar
  aicc <- NA
  if (n*p > npar + 2)
    aicc <- aic + (2*(npar+1)*(npar+2))/(n*p - npar - 2)
  bic  <- res.dev + npar*log(n*p)
  icl  <- bic
  outvect   <- this.fit$par
  alpha.hat <- outvect[1:(n-1)]
  alpha.hat <- c(alpha.hat,-sum(alpha.hat))
  beta.hat  <- outvect[n:(n+p-1)]
  k.hat     <- exp(outvect[(n+p):(n+2*p-1)])
  phi.hat   <- 1/k.hat
  lin.hat   <- outer(alpha.hat,beta.hat,FUN="+")
  mu.hat    <- exp(lin.hat)
  sigsq.hat <- mu.hat + mu.hat^2*matrix(phi.hat,n,p,byrow=T)
  k.hat[k.hat>1000] <- 1/0
  out1 <- round(c(n,p,max.ll,res.dev,npar,aic,aicc,bic,icl,max.ll),3)
  names(out1) <- c("n","p","Max.ll","Res.Dev.","npar","AIC","AICc",
                   "BIC","ICL","LLc")
  this.out <- list("info"=out1,
                   "alpha"=round(alpha.hat,3),
                   "beta"=round(beta.hat,3),
                   "phi"=round(phi.hat,3),
                   "k"=round(k.hat,3),
                   "mu"=round(mu.hat,3),
                   "sigma-sq"=round(sigsq.hat,3))
  this.out
}


# ---------------------------------------------------------------
#               Clustered models
# ---------------------------------------------------------------

# ---------------------------------------------------------------
#               Cluster Columns 
# ---------------------------------------------------------------

# Gamma (pattern) matrix is gamma(i,c).
# Row effects alpha(i) are subsumed into gamma matrix.

# Model C:
# --------

# Do not allow for column main effects beta(j).
# log mu(ijc) = gamma(ic) if col j in col group c.
# phi(j) = coeff. in quadratic for column j.
# sigsq(ij) = mu(ij) + phi(j) * mu(ij)^2
# Use k(j) = 1/phi(j), k is NB parameter.

# Independent parameters for optim:
# gamma.matrix, n by C  [1:(n*CG)]
# log(k)      [(n*CG+1):(n*CG+p)]

# Function to calculate -llc (excluding the kappa estimates).

NB.C.llc <- function(invect){
  gamma.mm <- matrix(invect[1:(n*CG)],n,CG)
  k.vv     <- exp(invect[(n*CG+1):(n*CG+p)])
  # Create linear array, n by p by CG:
  lin.aa <- array(NA,c(n,p,CG))
  for (jj in 1:p) lin.aa[,jj,] <- gamma.mm  
  # Create y, mu and k arrays:
  mu.aa <- exp(lin.aa)
  y.aa  <- array(y.mat,c(n,p,CG))
  k.aa  <- array(matrix(k.vv,n,p,byrow=T),c(n,p,CG))
  # Calculate array of complete log likelihood:
  llc.aa <- lgamma(y.aa+k.aa) - lgamma(y.aa+1) -
    lgamma(k.aa) + k.aa*log(k.aa) + y.aa*log(mu.aa) -
      (k.aa+y.aa)*log(k.aa+mu.aa)
  for (ii in 1:n) llc.aa[ii,,] <- llc.aa[ii,,]*ppc.in
  # Find llc:
  llc <- sum(llc.aa)
  llc
}


NB.C.fn <- function(CG=2,ppc.m=NULL,maxEMiter=10,maxreps=3){
  # If not given, set up ppc.m using kmeans:
  if (!is.matrix(ppc.m)){
    km <- kmeans(t(log(y.mat+1)),iter.max=50,CG,nstart=10)
    alloc <- km[[1]]
    ppc.m <- matrix(0,p,CG)
    for (jj in 1:p)
      ppc.m[jj,alloc[jj]] <- 1
  }
  # Set up starting parameters and do a preliminary M-step:
  kap.v   <- apply(ppc.m,2,mean)
  gamma.m <- matrix(0,n,CG)
  k.v     <- k0
  start.v <- c(gamma.m,log(k.v))
  CG      <<- CG
  ppc.in  <<- ppc.m
  this.fit <- optim(par = start.v,
                    fn = NB.C.llc,
                    method = "L-BFGS-B",
                    control = list(maxit=50,fnscale=-1))
  print("Initialising")
  print(paste("First optim value =",this.fit$value))
  # Update parameters:
  new.pars <- this.fit$par
  gamma.m  <- matrix(new.pars[1:(n*CG)],n,CG)
  k.v      <- exp(new.pars[(n*CG+1):(n*CG+p)])
  # Construct linear arrays, n by p by CG:
  lin.a <- array(NA,c(n,p,CG))
  for (jj in 1:p)
    lin.a[,jj,] <- gamma.m
  mu.a <- exp(lin.a)
  # Run the EM cycle until parameters have stabilised:
  print("Starting EM algorithm")
  EMoutvect <- c(gamma.m,k.v,kap.v)
  EMinvect  <- rep(1,length(EMoutvect))
  EMiter <- 1
  while(((EMiter==1)|(any(abs(EMinvect-EMoutvect)>1e-04)))&(EMiter<maxEMiter)){
    print(paste("EM iteration",EMiter))
    # E-step - Update posterior probabilities:
    for (jj in 1:p){
      # Do vector of log numerators (except kappa) for col jj:
      Ac.v <- rep(0,CG)
      for (cg in 1:CG){
        term <- sum(lgamma(y.mat[,jj]+k.v[jj])) -
          sum(lgamma(y.mat[,jj]+1)) -
            n*lgamma(k.v[jj]) +
              n*k.v[jj]*log(k.v[jj]) +
                sum(y.mat[,jj]*log(mu.a[,jj,cg])) -
                  sum((k.v[jj]+y.mat[,jj])*log(k.v[jj]+mu.a[,jj,cg]))
        Ac.v[cg] <- term
      }
      # Adjust to a maximum of zero:
      Bc.v <- Ac.v - max(Ac.v,na.rm=T)
      # Take exponential, multiply by kap.v:
      num.vect <- kap.v*exp(Bc.v)
      # Scale to add to 1:
      ppc.m[jj,] <- num.vect/sum(num.vect)
    }
    # M-step - Maximise LLC to update param estimates:
    kap.v   <- apply(ppc.m,2,mean)
    start.v <- c(gamma.m,log(k.v))
    ppc.in <<- ppc.m
    this.fit <- optim(par = start.v,
                      fn = NB.C.llc,
                      # ppc.in = ppc.m,
                      method = "L-BFGS-B",
                      control = list(maxit=50,fnscale=-1))
    # Optional extra runs:
    if (maxreps>0) for (it in 1:maxreps) if (this.fit$con != 0){
      start.v <- this.fit$par
      this.fit <- optim(par = start.v,
                        fn = NB.C.llc,
                        method = "L-BFGS-B",
                        hessian = T,
                        control = list(maxit=50,fnscale=-1))
  }
    print(paste("Current optim value =",round(this.fit$value,4)))
    # Update parameters:
    new.pars <- this.fit$par
    gamma.m  <- matrix(new.pars[1:(n*CG)],n,CG)
    k.v      <- exp(new.pars[(n*CG+1):(n*CG+p)])
    # Update arrays, n by p by CG:
    for (jj in 1:p)
      lin.a[,jj,] <- gamma.m
    mu.a <- exp(lin.a)
    # Update iteration number, vector of pars:
    EMiter <- EMiter + 1
    EMinvect <- EMoutvect
    EMoutvect <- c(gamma.m,k.v,kap.v)
  }
  # Find LLc:
  LLc <- this.fit$value + p*sum(xlogx.fn(kap.v))
  # Find overall log likelihood:
  entC <- -sum(xlogx.fn(ppc.m))
  logl <- LLc + entC
  # Save results:
  res.dev <- -2*logl
  npar <- length(start.v) + (CG-1) 
  aic  <- res.dev + 2*npar
  aicc <- NA
  if (n*p > npar + 2)
     aicc <- aic + (2*(npar+1)*(npar+2))/(n*p - npar - 2)
  bic <- -2*logl + npar*log(n*p)
  icl <- -2*LLc + npar*log(n*p)
  names(kap.v) <- paste("CG",1:CG,sep="")
  dimnames(gamma.m) <- list(rnames,paste("CG",1:CG,sep=""))
  dimnames(ppc.m) <- list(cnames,paste("CG",1:CG,sep=""))
  clus <- vector("list",CG)
  for (cg in 1:CG)
    clus[[cg]] <- cnames[ppc.m[,cg]==apply(ppc.m,1,max)]
  out1 <- round(c(n,p,logl,res.dev,npar,aic,aicc,bic,icl,LLc,n,CG),3)
  names(out1) <- c("n","p","Max.ll","Res.Dev.","npar","AIC","AICc",
                   "BIC","ICL","LLc","RG","CG")
  list("info"=out1,
       "kappa"=round(kap.v,3),
       "gamma"=round(gamma.m,3),
       "k"=round(k.v,3),
       "phi"=round(1/k.v,3),
       "ppc"=round(ppc.m,3),
       "column groups"=clus)
}

# -------------------------------------------------------------

# Model BC:
# ---------

# Allow for column main effects beta(j).
# log mu(ijc) = beta(j) + gamma(ic) if col j in col group c.
# phi(j) = coeff. in quadratic for column j.
# sigsq(ij) = mu(ij) + phi(j) * mu(ij)^2
# Use k(j) = 1/phi(j), k is NB parameter.

# Independent parameters for optim:
# beta vector [1:(p-1)]
# gamma.matrix, n by C  [(p):(p+n*CG-1)]
# log(k)      [(p+n*CG):(2*p+n*CG-1)]

# Function to calculate -llc (excluding the kappa estimates).

NB.BC.llc <- function(invect){
  beta.vv  <- invect[1:(p-1)]
  beta.vv <- c(beta.vv,-sum(beta.vv))
  gamma.mm <- matrix(invect[(p+0):(p+n*CG-1)],n,CG)
  k.vv     <- exp(invect[(p+n*CG+0):(2*p+n*CG-1)])
  # Create linear array, n by p by CG:
  lin.aa <- array(NA,c(n,p,CG))
  for (cc in 1:CG) lin.aa[,,cc] <-
    matrix(beta.vv,n,p,byrow=T) +
      matrix(gamma.mm[,cc],n,p)  
  # Create y, mu and k arrays:
  mu.aa <- exp(lin.aa)
  y.aa  <- array(y.mat,c(n,p,CG))
  k.aa  <- array(matrix(k.vv,n,p,byrow=T),c(n,p,CG))
  # Calculate array of complete log likelihood:
  llc.aa <- lgamma(y.aa+k.aa) - lgamma(y.aa+1) -
    lgamma(k.aa) + k.aa*log(k.aa) + y.aa*log(mu.aa) -
      (k.aa+y.aa)*log(k.aa+mu.aa)
  for (ii in 1:n) llc.aa[ii,,] <- llc.aa[ii,,]*ppc.in
  # Find llc:
  llc <- sum(llc.aa)
  llc
}


NB.BC.fn <- function(CG=2,ppc.m=NULL,maxEMiter=10,maxreps=3){
  # If not given, set up ppc.m using kmeans:
  if (!is.matrix(ppc.m)){
    km <- kmeans(t(log(y.mat+1)),iter.max=50,CG,nstart=10)
    alloc <- km[[1]]
    ppc.m <- matrix(0,p,CG)
    for (jj in 1:p)
      ppc.m[jj,alloc[jj]] <- 1
  }
  # Set up starting parameters:
  kap.v   <- apply(ppc.m,2,mean)
  gamma.m <- matrix(runif(n*CG),n,CG)
  beta.v  <- beta0
  k.v     <- k0
  # Do a preliminary M-step:
  start.v <- c(beta.v[-p],gamma.m,log(k.v))
  CG <<- CG
  print("Initialising:")
  ppc.in <<- ppc.m
  this.fit <- optim(par = start.v,
                    fn = NB.BC.llc,
                    # pp.in = ppc.m,
                    method = "L-BFGS-B",
                    control = list(maxit=50,fnscale=-1))
  print(paste("First optim value =",round(this.fit$value,4)))
  # Update parameters:
  new.pars <- this.fit$par
  beta.v   <- new.pars[1:(p-1)]
  beta.v <- c(beta.v,-sum(beta.v))
  gamma.m  <- matrix(new.pars[p:(p+n*CG-1)],n,CG)
  k.v      <- exp(new.pars[(p+n*CG):(2*p+n*CG-1)])
  # Construct linear arrays, n by p by CG:
  lin.a <- array(NA,c(n,p,CG))
  for (ii in 1:n) for (jj in 1:p) for (cg in 1:CG)
    lin.a[ii,jj,cg] <- beta.v[jj] + gamma.m[ii,cg]
  mu.a <- exp(lin.a)
  # Run the EM cycle until parameters have stabilised:
  print("Starting EM algorithm")
  EMoutvect <- c(beta.v,gamma.m,k.v,kap.v)
  EMinvect  <- rep(1,length(EMoutvect))
  EMiter <- 1
  while(((EMiter==1)|(any(abs(EMinvect-EMoutvect)>1e-04)))&(EMiter<maxEMiter)){
    print(paste("EM iteration",EMiter))
    # E-step - Update posterior probabilities:
    for (jj in 1:p){
      # Do vector of log numerators (except kappa) for col jj:
      Ac.v <- rep(0,CG)
      for (cg in 1:CG){
        term <- sum(lgamma(y.mat[,jj]+k.v[jj])) -
          sum(lgamma(y.mat[,jj]+1)) -
            n*lgamma(k.v[jj]) +
              n*k.v[jj]*log(k.v[jj]) +
                sum(y.mat[,jj]*log(mu.a[,jj,cg])) -
                  sum((k.v[jj]+y.mat[,jj])*log(k.v[jj]+mu.a[,jj,cg]))
        Ac.v[cg] <- term
      }
      # Adjust to a maximum of zero:
      Bc.v <- Ac.v - max(Ac.v,na.rm=T)
      # Take exponential, multiply by kap.v:
      num.vect <- kap.v*exp(Bc.v)
      # Scale to add to 1:
      ppc.m[jj,] <- num.vect/sum(num.vect)
    }
    # M-step - Maximise LLC to update param estimates:
    kap.v   <- apply(ppc.m,2,mean)
    start.v <- c(beta.v[-p],gamma.m,log(k.v))
    ppc.in <<- ppc.m
    this.fit <- optim(par = start.v,
                      fn = NB.BC.llc,
                      method = "L-BFGS-B",
                      control = list(maxit=50,fnscale=-1))
    # Optional extra runs:
    if (maxreps>0) for (it in 1:maxreps) if (this.fit$con != 0){
      start.v <- this.fit$par
      this.fit <- optim(par = start.v,
                        fn = NB.BC.llc,
                        # pp.in = pp.c,
                        method = "L-BFGS-B",
                        hessian = T,
                        control = list(maxit=50,fnscale=-1))
     }
    print(paste("Current optim value =",round(this.fit$value,4)))
    # Update parameters:
    new.pars <- this.fit$par
    beta.v   <- new.pars[1:(p-1)]
    beta.v <- c(beta.v,-sum(beta.v))
    gamma.m  <- matrix(new.pars[p:(p+n*CG-1)],n,CG)
    k.v      <- exp(new.pars[(p+n*CG):(2*p+n*CG-1)])
    # Update arrays, n by p by CG:
    for (ii in 1:n) for (jj in 1:p) for (cg in 1:CG)
      lin.a[ii,jj,cg] <- beta.v[jj] + gamma.m[ii,cg]
    mu.a <- exp(lin.a)
    # Update iteration number, vector of pars:
    EMiter <- EMiter + 1
    EMinvect <- EMoutvect
    EMoutvect <- c(beta.v,gamma.m,k.v,kap.v)
  }
  # Find LLc:
  LLc <- this.fit$value + p*sum(xlogx.fn(kap.v))
  # Find overall log likelihood:
  entC <- -sum(xlogx.fn(ppc.m))
  logl <- LLc + entC
  # Save results:
  res.dev <- -2*logl
  npar <- length(start.v) + (CG-1) 
  aic  <- res.dev + 2*npar
  aicc <- NA
  if (n*p > npar + 2)
     aicc <- aic + (2*(npar+1)*(npar+2))/(n*p - npar - 2)
  bic <- -2*logl + npar*log(n*p)
  icl <- -2*LLc + npar*log(n*p)
  names(kap.v) <- paste("CG",1:CG,sep="")
  dimnames(gamma.m) <- list(rnames,paste("CG",1:CG,sep=""))
  dimnames(ppc.m) <- list(cnames,paste("CG",1:CG,sep=""))
  clus <- vector("list",CG)
  for (cg in 1:CG)
    clus[[cg]] <- cnames[ppc.m[,cg]==apply(ppc.m,1,max)]
  out1 <- round(c(n,p,logl,res.dev,npar,aic,aicc,bic,icl,LLc,n,CG),3)
  names(out1) <- c("n","p","Max.ll","Res.Dev.","npar","AIC","AICc",
                   "BIC","ICL","LLc","RG","CG")
  list("info"=out1,
       "kappa"=round(kap.v,3),
       "beta"=round(beta.v,3),
       "gamma"=round(gamma.m,3),
       "k"=round(k.v,3),
       "phi"=round(1/k.v,3),
       "ppc"=round(ppc.m,3),
       "column groups"=clus)
}


# -------------------------------------------------------------

# Model ABC:
# ----------

# Allow for row and column main effects, alpha(i) and beta(j).
# log mu(ijc) = alpha(i) + beta(j) + gamma(ic) if col j in col group c.
# phi(j) = coeff. in quadratic for column j.
# sigsq(ij) = mu(ij) + phi(j) * mu(ij)^2
# Use k(j) = 1/phi(j), k is NB parameter.

# Same model as BC, but separate out alpha terms from the 
# gamma matrix, let delta matrix be the deviations.
# This makes the delta matrix show species composition patterns
# after allowing for different site means.
# Row sums of delta matrix are zero.

NB.ABC.fn <- function(BC.out){
  alpha.v <- apply(BC.out$gamma,1,mean)
  delta.m <- BC.out$gamma - alpha.v
  list("info"=BC.out$info,
       "kappa"=BC.out$kappa,
       "alpha"=round(alpha.v,3),
       "beta"=BC.out$beta,
       "delta"=round(delta.m,3),
       "k"=BC.out$k,
       "phi"=BC.out$phi,
       "ppc"=BC.out$ppc,
       "column groups"=BC.out$'column groups')
}


# -------------------------------------------------------------

