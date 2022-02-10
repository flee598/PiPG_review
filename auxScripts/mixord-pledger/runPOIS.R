# runPOIS.R

# Count data

# Code to fit  models to count data.
# This code uses functions from   funPOIS.R
# Run the code interactively by cutting and pasting blocks of
# commands into R.

# Need to define y.mat,  the input matrix of count data.

# -----------------------------------------------------------------
#             Define the data matrix
# -----------------------------------------------------------------

# The data should be a matrix of count data, with each site
# being a row and each species being a column.
# It should be called   y.mat
# We use the spider data as an example.

source("spiderin.R")

y.mat <- spider.mat            # 28 by 12 matrix of counts        

# ------------------------------------------------------
#                   Setting up 
# ------------------------------------------------------
 
n       <- nrow(y.mat)           # Number of sites
p       <- ncol(y.mat)           # Number of species
lambda0 <- mean(y.mat)           # Overall mean
FC      <- sum(lgamma(y.mat+1))  # Factorial constant
rnames  <- rownames(y.mat)
cnames  <- colnames(y.mat)

# Read in the functions:

source("funPOIS.R")



# ------------------------------------------------------
#             Null Model 
# ------------------------------------------------------

# mu = E(Y)
# log(mu(i,j))= theta

POIS.Null.out <- POIS.Null.fn()

print(POIS.Null.out)


# ------------------------------------------------------
#             Model A 
# ------------------------------------------------------

# Allow for site effects.
# Species effects not allowed for in mu.
# mu = E(Y)
# log(mu(i,j))= alpha(i)

POIS.A.out <- POIS.A.fn()

print(POIS.A.out)



# ------------------------------------------------------
#             Model B 
# ------------------------------------------------------

# Allow for species effects.
# Site effects not allowed for in mu.
# mu = E(Y)
# log(mu(i,j))= beta(j)

POIS.B.out <- POIS.B.fn()

print(POIS.B.out)



# ------------------------------------------------------
#             Model AB 
# ------------------------------------------------------

# Adjust for both site and species effects.
# No-association model.
# mu = E(Y)
# log(mu(i,j))= alpha(i) + beta(j)  with sum(alpha) = 0

POIS.AB.out <- POIS.AB.fn()

print(POIS.AB.out)



# ------------------------------------------------------
#     Model C  Rows separate, column clustered
# ------------------------------------------------------

# Ordination from Model C is driven by species commonness,
# site richness and species composition (turnover).
# log(lambda(ijc) = gamma(ic)
# The gamma matrix provides coordinates for the ordination.


# Set up matrix of residuals for a starting point inside POIS.C.fn

res.mat <- y.mat/POIS.AB.out$mu

# Two column groups (set by "CG=2"):

POIS.C2.out <- POIS.C.fn(CG=2)
   
print(POIS.C2.out)

# Plot ordination of sites:

plot(POIS.C2.out$lambda,type="n",
     xlab="Dimension.1",
     ylab="Dimension.2",
     main="Spiders, ordination of sites, POIS, model C2")
text(POIS.C2.out$lambda,as.character(1:n))


# For model comparison using AIC, can do models with more column groups:

POIS.C3.out <- POIS.C.fn(CG=3)   
print(POIS.C3.out)

POIS.C4.out <- POIS.C.fn(CG=4)   
print(POIS.C4.out)

POIS.C5.out <- POIS.C.fn(CG=5)   
print(POIS.C5.out)



# ------------------------------------------------------
#             Model BC 
# ------------------------------------------------------

# Model BC allows for common versus rare species, so ordination
# from Model BC is determined by rich versus poor sites and by
# species turnover.
# Allow for species effects.
# log(lambda(ijc) = beta(j) + gamma(ic), with sum(beta) = 0

beta0 <- POIS.AB.out$beta
beta0 <- beta0 - mean(beta0)

# Two column groups, random starts:

POIS.BC2.out <- POIS.BC.fn(CG=2)
   
print(POIS.BC2.out)

# Plot ordination of rows:

plot(POIS.BC2.out$gamma,type="n",
     xlab="Dimension.1",
     ylab="Dimension.2",
     main="Spiders, ordination of sites, POIS, model BC2")
text(POIS.BC2.out$gamma,as.character(1:n))


# For model comparison using AIC, can do models with more column groups:

POIS.BC3.out <- POIS.BC.fn(CG=3)   
print(POIS.BC3.out)

POIS.BC4.out <- POIS.BC.fn(CG=4)   
print(POIS.BC4.out)

POIS.BC5.out <- POIS.BC.fn(CG=5)   
print(POIS.BC5.out)


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

POIS.ABC3.out <- POIS.ABC.fn(POIS.BC3.out)
print(POIS.ABC3.out)


# The points (rows) in POIS.ABC3.out$delta are coplanar,
# all on the plane x+y+z=0.
# Need to make this a 2D plot.

x <- POIS.ABC3.out$delta-min(POIS.ABC3.out$delta)
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
     main="Spiders, ordination of sites, POIS, model ABC3")
#     sub="Site:  site effect separated out, 
#          ordination by species composition")
text(points2D,as.character(1:n))











