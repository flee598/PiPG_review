# funPOIS.R 

# Functions for count data, pattern detection models.
# Poisson building block.
# Column clustering.

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

# GLM with log mu(ij) = theta (constant).
# mu = lambda

POIS.Null.fn <- function(){
  # Estimation:
  lam <- mean(y.mat)
  thet <- log(lam)
  logl <- - sum(y.mat) + sum(y.mat)*thet - FC
  # Save results:
  res.dev <- -2*logl
  npar <- 1
  aic  <- res.dev + 2*npar
  aicc    <- aic + (2*(npar+1)*(npar+2))/(n*p - npar - 2)
  bic     <- -2*logl + npar*log(n*p)
  out1    <-  round(c(n,p,logl,res.dev,npar,aic,aicc,bic,bic,logl,1,1),3)
  names(out1) <- c("n","p","Max.ll","Res.Dev","npar","AIC","AICc","BIC",
                   "ICL","LLc","RG","CG")
  list("info"=out1,"lambda"=round(lam,3),theta=round(thet,3))
}


# ------------ Model A ----------------------------------------------

# GLM with log mu(ij) = alpha(i).
# Linear predictor: log(mu(ij)) = alpha(i)

POIS.A.fn <- function(){
  # Estimation:
  lambda.v <- apply(y.mat,1,mean)
  alpha.v  <- log(lambda.v)
  logl <- - sum(y.mat) + sum(apply(y.mat,1,sum)*alpha.v) - FC
  # Save results:
  res.dev <- -2*logl
  npar    <- n
  aic     <- res.dev + 2*npar
  aicc    <- aic + (2*(npar+1)*(npar+2))/(n*p - npar - 2)
  bic     <- -2*logl + npar*log(n*p)
  out1    <-  round(c(n,p,logl,res.dev,npar,aic,aicc,bic,bic,logl,n,1),3)
  names(out1) <- c("n","p","Max.ll","Res.Dev","npar","AIC","AICc","BIC",
                   "ICL","LLc","RG","CG")
  list("info"=out1,
       "lambda"=round(lambda.v,3),
       "alpha"=round(alpha.v,3))
}


# ------------ Model B --------------------------------------------

# GLM with log mu(ij) = beta(j).
# Linear predictor: log(mu(ij)) = beta(j)

POIS.B.fn <- function(){
  # Estimation:
  lambda.v <- apply(y.mat,2,mean)
  beta.v   <- log(lambda.v)
  logl <- - sum(y.mat) + sum(apply(y.mat,2,sum)*beta.v) - FC
  # Save results:
  res.dev <- -2*logl
  npar <- p
  aic  <- res.dev + 2*npar
  aicc    <- aic + (2*(npar+1)*(npar+2))/(n*p - npar - 2)
  bic     <- -2*logl + npar*log(n*p)
  out1    <-  round(c(n,p,logl,res.dev,npar,aic,aicc,bic,bic,logl,1,p),3)
  names(out1) <- c("n","p","Max.ll","Res.Dev","npar","AIC","AICc","BIC",
                   "ICL","LLc","RG","CG")
  list("info"=out1,
       "lambda"=round(lambda.v,3),
       "beta"=round(beta.v,3))
}



# ------------ Model AB -----------------------------------------

# Basic null model of no association.
# log mu(ij) = alpha(i) + beta(j) with constraint sum alpha(i)=0.

# Parameters for optim:
# alpha[-n]   [1:(n-1)]
# beta        [n:(n+p-1)]

POIS.AB.ll <- function(invect){
  # Read in parameters:
  alpha.v <- invect[1:(n-1)]
  alpha.v <- c(alpha.v,0-sum(alpha.v))
  beta.v  <- invect[n:(n+p-1)]
  # Find the linear predictor:
  lin.m   <- outer(alpha.v,beta.v,FUN="+")
  # Find mu:
  mu.m     <- exp(lin.m)
  # Calculate NB likelihood:
  ll <- -1000000
  if (min(c(mu.m))>0) ll <-
    sum(-mu.m + y.mat*log(mu.m)) - FC
  ll
}

POIS.AB.fn <- function(start.v=rep(0,n+p-1),maxreps=3){
  # Run the optimisation, with repeats if necessary:
  this.fit <- optim(par=start.v,
                    fn=POIS.AB.ll,
                    method = "L-BFGS-B",
                    control = list(fnscale=-1))
  if (maxreps>0) 
  for (i in 1:maxreps) if (this.fit$con != 0){
    start.vect <- this.fit$par
    this.fit <- optim(par=start.vect,
                      fn=POIS.AB.ll,
                      method = "L-BFGS-B",
                    control = list(fnscale=-1))
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
  lin.hat   <- outer(alpha.hat,beta.hat,FUN="+")
  mu.hat    <- exp(lin.hat)
  out1 <- round(c(n,p,max.ll,res.dev,npar,aic,aicc,bic,icl,max.ll),3)
  names(out1) <- c("n","p","Max.ll","Res.Dev.","npar","AIC","AICc",
                   "BIC","ICL","LLc")
  this.out <- list("info"=out1,
                   "alpha"=round(alpha.hat,3),
                   "beta"=round(beta.hat,3),
                   "mu"=round(mu.hat,3))
  this.out
}




# ------------ Model C ----------------------------------------

# log(lambda(ijc)) = gamma(ic).
# Allow for large counts.

POIS.C.fn <- function(CG,maxEMiter=20,kstarts=20,start.ppc=0){
  if (is.matrix(start.ppc)){
    ppc.m  <- start.ppc
    } else{
    # Generate start from kmeans on y.mat:
    temp <- kmeans(t(y.mat),CG,nstart=kstarts)$cluster
    ppc.m <- matrix(0,p,CG)
    for (cg in 1:CG) for (j in 1:p) 
      if (temp[j]==cg) ppc.m[j,cg] <- 1
    }
  # Find starting parameters:
  kap.v <- apply(ppc.m,2,mean)
  lam.m <- y.mat %*% ppc.m
  for (cg in 1:CG) 
    lam.m[,cg] <- lam.m[,cg]/p/kap.v[cg]
  # Calculate lambda array, n by p by CG:
  lam.arr <- array(NA,c(n,p,CG))
  for (jj in 1:p)
    lam.arr[,jj,] <- lam.m 
  # Run the EM cycle until parameters have stabilised:
  print("Starting EM algorithm")
  EMoutvect <- c(lam.m,kap.v)
  EMinvect  <- rep(1,length(EMoutvect))
  EMiter <- 1
  while(((EMiter==1)|(any(abs(EMinvect-EMoutvect)>1e-04)))&(EMiter<maxEMiter)){
    print(paste("EM iteration",EMiter))
    # E-step - Update posterior probabilities:
    for (jj in 1:p){
      # Do log numerators (except kappa) for col jj:
      Ac.v <- rep(0,CG)
      for (cg in 1:CG) if (kap.v[cg]>0){
        term <- - sum(lam.arr[,jj,cg])
        for (ii in 1:n) if (lam.arr[ii,jj,cg]>0)
          term <- term + 
                  y.mat[ii,jj]*log(lam.arr[ii,jj,cg])
          Ac.v[cg] <- term
      }
      # Adjust to a maximum of zero:
      Bc.v <- Ac.v - max(Ac.v,na.rm=T)
      # Take exponential, multiply by kap.v:
      num.vect <- kap.v*exp(Bc.v)
      # Scale to add to 1:
      ppc.m[jj,] <- num.vect/sum(num.vect)
    }
    kap.v  <- apply(ppc.m,2,mean)
    # M-step - Update parameter estimates:
    kap.v <- apply(ppc.m,2,mean)
    lam.m <- y.mat %*% ppc.m
    for (cg in 1:CG) 
      lam.m[,cg] <- lam.m[,cg]/p/kap.v[cg]
    # Calculate lambda array, n by p by CG:
    lam.arr <- array(NA,c(n,p,CG))
    for (jj in 1:p)
      lam.arr[,jj,] <- lam.m 
    EMiter    <- EMiter + 1
    EMinvect  <- EMoutvect
    EMoutvect <- c(lam.m,kap.v)
  }
  # Find LLc:
  LLc <- p*sum(xlogx.fn(kap.v)) - FC
  for (ii in 1:n) for (jj in 1:p) for (cg in 1:CG)
    if (lam.arr[ii,jj,cg]>0) LLc <- LLc +
      ppc.m[jj,cg]*(y.mat[ii,jj]*log(lam.arr[ii,jj,cg]) -
                   lam.arr[ii,jj,cg])
  # Calculate log likelihood using entropy:
  entC <- -sum(xlogx.fn(ppc.m))
  logl <- LLc + entC
  # Calculate gamma matrix = log(lambda matrix):
  gam.m <- log(lam.m)
  # Save results:
  res.dev <- -2*logl
  npar <- CG*n + CG - 1  
  aic  <- res.dev + 2*npar
  aicc <- NA
  if (n*p > npar + 2)
     aicc <- aic + (2*(npar+1)*(npar+2))/(n*p - npar - 2)
  bic <- -2*logl + npar*log(n*p)
  icl <- -2*LLc + npar*log(n*p)
  names(kap.v) <- paste("CG",1:CG,sep="")
  dimnames(gam.m) <- list(rnames,paste("CG",1:CG,sep=""))
  dimnames(ppc.m) <- list(cnames,paste("CG",1:CG,sep=""))
  clus <- vector("list",CG)
  for (cg in 1:CG)
    clus[[cg]] <- cnames[ppc.m[,cg]==apply(ppc.m,1,max)]
  out1 <- round(c(n,p,logl,res.dev,npar,aic,aicc,bic,icl,LLc,n,CG),3)
  names(out1) <- c("n","p","Max.ll","Res.Dev.","npar","AIC","AICc",
                   "BIC","ICL","LLc","RG","CG")
  list("info"=out1,
       "kappa"=round(kap.v,3),
       "lambda"=round(lam.m,3),
       "gamma"=round(gam.m,3),
       "ppc"=round(ppc.m,3),
       "column groups"=clus)
}


# ------------ Model BC ----------------------------------------

# Adjust for column sums only, then do column clusters.
# log(lambda(ijc)) = beta(j) + gamma(ic)
# Sum  beta(j) = 0
# Have all rows varying (n levels).
# Allow for large counts.
# Vector of pars, beta[1:(p-1)], gamma[p:(p-1*n*CG)]

POIS.BC.llc <- function(invect){
  # This function omits the kappa term.
  beta.vv  <- invect[1:(p-1)]
  beta.vv <- c(beta.vv,-sum(beta.vv))
  gamma.mm <- matrix(invect[(p+0):(p+n*CG-1)],n,CG)
  # Create linear array, n by p by CG:
  lin.aa <- array(NA,c(n,p,CG))
  for (cc in 1:CG) lin.aa[,,cc] <-
    matrix(beta.vv,n,p,byrow=T) +
      matrix(gamma.mm[,cc],n,p)  
  # Create y and lambda arrays:
  lam.aa <- exp(lin.aa)
  y.aa  <- array(y.mat,c(n,p,CG))
  # Calculate array of complete log likelihood:
  llc.aa <- - lam.aa + y.aa*log(lam.aa)
  for (ii in 1:n) llc.aa[ii,,] <- llc.aa[ii,,]*ppc.in
  # Find llc:
  llc <- sum(llc.aa) - FC
  llc
}


POIS.BC.fn <- function(CG=2,ppc.m=NULL,maxEMiter=10,maxreps=3){
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
  # Do a preliminary M-step:
  start.v <- c(beta.v[-p],gamma.m)
  CG <<- CG
  print("Initialising:")
  ppc.in <<- ppc.m
  this.fit <- optim(par = start.v,
                    fn = POIS.BC.llc,
                    method = "L-BFGS-B",
                    control = list(maxit=50,fnscale=-1))
  print(paste("First optim value =",round(this.fit$value,4)))
  # Update parameters:
  new.pars <- this.fit$par
  beta.v   <- new.pars[1:(p-1)]
  beta.v <- c(beta.v,-sum(beta.v))
  gamma.m  <- matrix(new.pars[p:(p+n*CG-1)],n,CG)
  # Construct linear arrays, n by p by CG:
  lin.a <- array(NA,c(n,p,CG))
  for (ii in 1:n) for (jj in 1:p) for (cg in 1:CG)
    lin.a[ii,jj,cg] <- beta.v[jj] + gamma.m[ii,cg]
  lam.a <- exp(lin.a)
  # Run the EM cycle until parameters have stabilised:
  print("Starting EM algorithm")
  EMoutvect <- c(beta.v,gamma.m,kap.v)
  EMinvect  <- rep(1,length(EMoutvect))
  EMiter <- 1
  while(((EMiter==1)|(any(abs(EMinvect-EMoutvect)>1e-04)))&(EMiter<maxEMiter)){
    print(paste("EM iteration",EMiter))
    # E-step - Update posterior probabilities:
    for (jj in 1:p){
      # Do vector of log numerators (except kappa and log(y(ij)!) for col jj:
      Ac.v <- rep(0,CG)
      for (cg in 1:CG){
        term <- sum(-lam.a[,jj,cg]) + sum(y.mat[,jj]*log(lam.a[,jj,cg]))
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
    start.v <- c(beta.v[-p],gamma.m)
    ppc.in <<- ppc.m
    this.fit <- optim(par = start.v,
                      fn = POIS.BC.llc,
                      method = "L-BFGS-B",
                      control = list(maxit=50,fnscale=-1))
    # Optional extra runs:
    if (maxreps>0) for (it in 1:maxreps) if (this.fit$con != 0){
      start.v <- this.fit$par
      this.fit <- optim(par = start.v,
                        fn = POIS.BC.llc,
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
    # Update arrays, n by p by CG:
    for (ii in 1:n) for (jj in 1:p) for (cg in 1:CG)
      lin.a[ii,jj,cg] <- beta.v[jj] + gamma.m[ii,cg]
    lam.a <- exp(lin.a)
    # Update iteration number, vector of pars:
    EMiter <- EMiter + 1
    EMinvect <- EMoutvect
    EMoutvect <- c(beta.v,gamma.m,kap.v)
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
       "ppc"=round(ppc.m,3),
       "column groups"=clus)
}



# ------------ Model ABC ---------------------------------

# log lambda(ijc) = alpha(i) + beta(j) + gamma(ic) if col j in col group c.

# Same model as BC, but separate out alpha terms from the 
# gamma matrix, let delta matrix be the deviations.
# Partition gamma(ic) from Model BC into alpha(i) and delta(ic).
# alpha(i) are row sums from gamma(ic).
# delta(ic) are deviations from the row sums.
# This makes the delta matrix show species composition patterns
# after allowing for different site means.
# Row sums of delta matrix are zero.

POIS.ABC.fn <- function(POIS.BC.out){
  alpha.v <- apply(POIS.BC.out$gamma,1,mean)
  delta.m <- POIS.BC.out$gamma - alpha.v
  list("info"=POIS.BC.out$info,
       "kappa"=POIS.BC.out$kappa,
       "alpha"=round(alpha.v,3),
       "beta"=POIS.BC.out$beta,
       "delta"=round(delta.m,3),
       "ppc"=POIS.BC.out$ppc,
       "column groups"=POIS.BC.out$'column groups')
}


# -------------------------------------------------------------


