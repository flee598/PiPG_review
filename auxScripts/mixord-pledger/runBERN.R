# runBERN.R

# Binary data

# Code to fit models to binary data.
# This code uses functions from   funBERN.R
# Run the code interactively by cutting and pasting blocks of
# commands into R.

# Need to define y.mat,  the input matrix of binary data.

# -----------------------------------------------------------------
#             Define the data matrix
# -----------------------------------------------------------------

# The data should be a matrix of binary data, with each site
# being a row and each species being a column.
# It should be called   y.mat
# We use the Tikus coral data as an example.

Tikus.df <- read.csv("TikusData.csv",T)
Tikus.mat <- as.matrix(Tikus.df[,-1])
rownames(Tikus.mat) <- Tikus.df[,1]

y.mat <- Tikus.mat            # 19 by 18 matrix of binary data        

# ------------------------------------------------------
#                   Setting up 
# ------------------------------------------------------
 
n       <- nrow(y.mat)           # Number of sites
p       <- ncol(y.mat)           # Number of species
theta0  <- mean(y.mat)           # Overall mean
rnames  <- rownames(y.mat)
cnames  <- colnames(y.mat)
rs.vect <- apply(y.mat,1,sum)    # Row sums
cs.vect <- apply(y.mat,2,sum)    # Column sums

# Read in the functions:

source("funBERN.R")



# ------------------------------------------------------
#             Null Model 
# ------------------------------------------------------

# mu = E(Y)
# logit(mu(i,j))= theta

BERN.Null.out <- BERN.Null.fn()

BERN.Null.out


# ------------------------------------------------------
#             Model A 
# ------------------------------------------------------

# Allow for site effects.
# Species effects not allowed for in mu.
# mu = E(Y)
# logit(mu(i,j))= alpha(i)

BERN.A.out <- BERN.A.fn()

BERN.A.out


# ------------------------------------------------------
#             Model B 
# ------------------------------------------------------

# Allow for species effects.
# Site effects not allowed for in mu.
# mu = E(Y)
# logit(mu(i,j))= beta(j)

BERN.B.out <- BERN.B.fn()

BERN.B.out


# ------------------------------------------------------
#             Model AB 
# ------------------------------------------------------

# Adjust for both site and species effects.
# No-association model.
# mu = E(Y)
# logit(mu(i,j))= alpha(i) + beta(j)  with sum(alpha) = 0

BERN.AB.out <- BERN.AB.fn()

BERN.AB.out



# ------------------------------------------------------
#     Model C  Rows separate, column clustered
# ------------------------------------------------------

# Ordination from Model C is driven by species commonness,
# site richness and species composition (turnover).
# log(lambda(ijc) = gamma(ic)
# The gamma matrix provides coordinates for the ordination.

# Two column groups (set by "CG=2").

BERN.C2.out <- BERN.C.fn(CG=2)
   
BERN.C2.out


# Plot ordination of sites:

plot(BERN.C2.out$mu,type="n",
     xlab="Dimension.1",
     ylab="Dimension.2",
     main="Corals, ordination of sites, BERN, model C2")
text(BERN.C2.out$mu,as.character(1:n))

# For model comparison using AIC, can do models with more column groups:

BERN.C3.out <- BERN.C.fn(CG=3)

BERN.C3.out


BERN.C4.out <- BERN.C.fn(CG=4)

BERN.C4.out


BERN.C5.out <- BERN.C.fn(CG=5)

BERN.C5.out



# ------------------------------------------------------
#             Model BC 
# ------------------------------------------------------

# Model BC allows for common versus rare species, so ordination
# from Model BC is determined by rich versus poor sites and by
# species turnover.
# Allow for species effects.
# log(lambda(ijc) = beta(j) + gamma(ic), with sum(beta) = 0

beta0 <- BERN.AB.out$beta
beta0 <- beta0 - mean(beta0)

# Two column groups, random starts.

BERN.BC2.out <- BERN.BC.fn(CG=2)
   
BERN.BC2.out


# Plot ordination of rows:

plot(BERN.BC2.out$gamma,type="n",
     xlab="Dimension.1",
     ylab="Dimension.2",
     main="Corals, ordination of sites, BERN, model BC2")
text(BERN.BC2.out$gamma,as.character(1:n))

# For model comparison using AIC, can do models with more column groups:

BERN.BC3.out <- BERN.BC.fn(CG=3)   

BERN.BC3.out


BERN.BC4.out <- BERN.BC.fn(CG=4)   

BERN.BC4.out


BERN.BC5.out <- BERN.BC.fn(CG=5)   

BERN.BC5.out


# ------------------------------------------------------
#             Model ABC 
# ------------------------------------------------------

# log lambda(ijc) = alpha(i) + beta(j) + gamma(ic) if col j in col group c.

# Model ABC allows for rich versus poor sites and common versus rare
# species, so ordination from Model ABC is determined solely by species
# composition.

# Same model as BC, but separate out alpha terms from the 
# gamma matrix, let delta matrix be the deviations.
# Partition gamma(ic) from Model BC into alpha(i) and delta(ic).
# alpha(i) are row sums from gamma(ic).
# delta(ic) are deviations from the row sums.
# This makes the delta matrix show species composition patterns
# after allowing for different site means.
# Row sums of delta matrix are zero.

# Model with 3 column groups, in order to obtain 2D plot:

BERN.ABC3.out <- BERN.ABC.fn(BERN.BC3.out)

BERN.ABC3.out


# The points (rows) in NB.ABC3.out$gamma are coplanar,
# all on the plane x+y+z=0.
# Need to make this a 2D plot.

x <- BERN.ABC3.out$delta-min(BERN.ABC3.out$delta)
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
     main="Corals, ordination of sites, BERN, model ABC3")
text(points2D,as.character(1:n))


