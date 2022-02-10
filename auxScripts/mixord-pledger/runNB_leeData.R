# runNB.R

# Count data 
# NB = negative binomial

# This code uses functions from   funNB.R.
# Run the code interactively by cutting and pasting blocks of
# commands into R.

# -----------------------------------------------------------------
#             Define the data matrix
# -----------------------------------------------------------------

# The data should be a matrix of count data, with each site
# being a row and each species being a column.
# It should be called   y.mat
# We use the spider data as an example.

source("spiderin.R")

y.mat <- spider.mat            # Matrix of counts        

inverts.mat <- as.matrix(inverts)
y.mat <- inverts.mat

# 28 samples (rows), 12 species (columns).

sitenames <- rownames(inverts)
spnames   <- colnames(inverts.mat)
spnam     <- colnames(inverts.mat)

# ------------------------------------------------------
#                   Setting up 
# ------------------------------------------------------
 
n       <- nrow(y.mat)           # Number of sites
p       <- ncol(y.mat)           # Number of species
lambda0 <- mean(y.mat)           # Overall mean.
rnames  <- rownames(y.mat)
cnames  <- colnames(y.mat)

# Read in the functions:

source("funNB.R")


# ------------------------------------------------------
#             Null Model 
# ------------------------------------------------------

# mu = E(Y)
# log(mu(i,j))= lambda

NB.Null.out <- NB.Null.fn()

NB.Null.out


# ------------------------------------------------------
#             Model A 
# ------------------------------------------------------

# Allow for site effects.
# Species effects not allowed for in mu.
# mu = E(Y)
# log(mu(i,j))= alpha(i)

NB.A.out <- NB.A.fn()

NB.A.out

# Species 11, T.ter, has very high k.hat. 
# It's the only k above 1.
# This column has the lowest variance; trying to make var = mean. 
# Other estimates and info output looks OK.



# ------------------------------------------------------
#             Model B 
# ------------------------------------------------------

# Allow for species effects.
# Site effects not allowed for in mu.
# mu = E(Y)
# log(E[Y(i,j)]) = log(mu(i,j))= beta(j)

NB.B.out <- NB.B.fn()

NB.B.out



# ------------------------------------------------------
#             Model AB 
# ------------------------------------------------------

# Adjust for both site and species effects.
# No-association model.
# mu = E(Y)
# log(E[Y(i,j)]) = log(mu(i,j))= alpha(i) + beta(j)  with sum(alpha) = 0

NB.AB.out <- NB.AB.fn()

NB.AB.out

# Same concern as with model A, 11th k large.
# Trying to make variance = mean.


# Save some pars for later starts:

temp.out <- NB.AB.out

alpha0 <- temp.out$alpha
beta0  <- temp.out$beta
k0     <- temp.out$k
k0[k0>1000] <- 1000



# ------------------------------------------------------
#             Model C 
# ------------------------------------------------------

# Ordination from Model C is driven by species commonness,
# site richness and species composition (turnover).
# Do not allow for separate row or column main effects (alpha, beta).
# Row  main effects (alpha) are subsumed into the pattern matrix, gamma.
# log(mu(ij|c))= gamma(ic).
# phi(j) = coeff. in quadratic for column j.
# sigsq(ij) = mu(ij) + mu(ij)^2*phi(j).
# Use k(j) = 1/phi(j).

# Independent parameters for optim:
# gamma.matrix, n by C  [1:(n*CC)]
# log(k)      [(n*CC+1):(n*CC+p)]

# Default is CG=2,maxEMiter=10, maxreps=3.
# Two column groups (set by "CG=2"):

NB.C2.out <- NB.C.fn(CG=2)

NB.C2.out


# Plot ordination of sites:

plot(NB.C2.out$gamma,type="n",
     xlab="Dimension.1",
     ylab="Dimension.2",
     main="Spiders, ordination of sites, Neg.Bin., model C2")
text(NB.C2.out$gamma,as.character(1:n))


# For model comparison using AIC, can do models with more column groups:

NB.C3.out <- NB.C.fn(CG=3)

NB.C3.out


NB.C4.out <- NB.C.fn(CG=4)

NB.C4.out


NB.C5.out <- NB.C.fn(CG=5)

NB.C5.out



# ------------------------------------------------------
#             Model BC 
# ------------------------------------------------------

# Model BC allows for common versus rare species, so ordination
# from Model BC is determined by rich versus poor sites and by
# species turnover.
# Site main effects are subsumed into the pattern matrix, gamma.
# mu(i,j,c) = E(Y(i,j)) given col j in col group c.
# log(E[Y(i,j)] | j in c) = log(mu(i,j,c))= beta(j) + gamma(ic).
# phi(j) = coeff. in quadratic for column j.
# sigsq(ij) = mu(ij) + mu(ij)^2*phi(j)

# Independent parameters for optim:
# beta vector [1:p]
# gamma.matrix, n by C  [(p+1):(p+n*CC)]
# log(k)      [(p+n*CC+1):(2*p+n*CC)]

# Default is CG=2,maxEMiter=10, maxreps=3.
# Two column groups:

NB.BC2.out <- NB.BC.fn(CG=2)

NB.BC2.out


# Plot ordination of sites:

plot(NB.BC2.out$gamma,type="n",
     xlab="Dimension 1",
     ylab="Dimension 2",
     main="Spiders, ordination of sites, Neg.Bin., model BC2")
text(NB.BC2.out$gamma,as.character(1:n))


# For model comparison using AIC, can do models with more column groups:

NB.BC3.out <- NB.BC.fn(CG=3)

NB.BC3.out


NB.BC4.out <- NB.BC.fn(CG=4)

NB.BC4.out


NB.BC5.out <- NB.BC.fn(CG=5)

NB.BC5.out



# ------------------------------------------------------
#     Model ABC - partition matrix from BC
# ------------------------------------------------------

# Model ABC allows for rich versus poor sites and common versus rare
# species, so ordination from Model ABC is determined solely by species
# composition.

# Same model as BC, but separate out alpha terms from the gamma matrix.
# This makes the gamma matrix show species composition patterns
# after allowing for different site means.

# Model with 3 column groups, in order to obtain 2D plot:

NB.ABC3.out <- NB.ABC.fn(NB.BC3.out)
print(NB.ABC3.out)


# The points (rows) in NB.ABC3.out$gamma are coplanar,
# all on the plane x+y+z=0.
# Need to make this a 2D plot.

x <- NB.ABC3.out$delta-min(NB.ABC3.out$delta)
s <- rowSums(x)
x <- x/s
top <- sqrt(3)/2
xlim <- c(-0.03, 1.03)
ylim <- c(-1, top)
xp <- x[, 2] + x[, 3]/2
yp <- x[, 3] * top
points2D <- cbind(xp,yp)

plot(points2D,type="n",
     xlab="Dimension 1",
     ylab="Dimension 2",
     main="Spiders, ordination of sites, Neg.Bin, model ABC3")
#     sub="Site:  site effect separated out, 
#          ordination by species composition")
text(points2D,as.character(1:n))

