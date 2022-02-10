# See: Hui, F.K.C. (2016). boral - Bayesian Ordination and Regression Analysis of Multivariate Abundance Data in r. Methods in Ecology and Evolution, 7, 744–750.


library(boral)
library(gllvm)
library(tidyverse)
data(spider)



inverts <- read_csv(file = 'Data/inverts.csv')
chem <- read_csv(file = 'Data/physicochemical.csv')
locations <- read_csv(file = '../../Data/Site_Coords_WGS84.csv')


y <- spider$abun
y <- inverts[,-1]

fit.lvmp <- boral(y = y, family = "poisson", lv.control = list(num.lv = 2), row.eff = "fixed") # Poisson model

summary(fit.lvmp)
plot(fit.lvmp, ask = F, mfrow = c(2,2))
par(mfrow = c(1,2))
lvsplot(fit.lvmp) # Biplot of latent variable medians from Poisson model
lvsplot(fit.lvmp, ind.spp = 5)


###  GLLVM
y <- inverts[,-1]
X <- chem %>%
  filter(site %in% inverts$Site)

X <- scale(X[,-c(1, 35)])  # scale and drop site and landform

fit.env.p <- gllvm(y, X, formula = ~ ProportionSilt + ProportionMacrophyte + Elevation + DO + ProportionRiparianCover0.5, family = poisson())


fit.env.nb <- gllvm(y, X, formula = ~ ProportionSilt + ProportionMacrophyte + Elevation + DO + ProportionRiparianCover0.5, family = 'negative.binomial') #AICc =  4445.57

fit.env.nb2 <- gllvm(y, X, formula = ~ Pasture, family = 'negative.binomial') #AICc =  5137.693


ordiplot.gllvm(fit.env.nb, biplot = TRUE, ind.spp = 10)
ordiplot.gllvm(fit.env.nb2, biplot = TRUE, ind.spp = 10)

coefplot(fit.env.nb, which.Xcoef = 1:2, cex.ylab = 0.7, xlim.list = list(NULL, NULL, c(-4, 4)))

par(mfrow = c(1,2))
ordiplot.gllvm(fit.env.nb2, biplot = TRUE, ind.spp = 10)
coefplot(fit.env.nb2, which.Xcoef = 1, cex.ylab = 0.7, xlim.list = list(NULL, NULL, c(-4, 4)))

criterias <-NULL 
for(i in 1:5) {
  fiti <- gllvm(y, X,family ="negative.binomial",num.lv =i,sd.errors =FALSE, formula = ~ .,seed =1234) 
  criterias[i + 1] <- summary(fiti)$AICc 
  names(criterias)[i + 1] =i }



cr <- getResidualCor(fit.env.nb2)
library(corrplot)

corrplot(cr[order.single(cr), order.single(cr)], diag = FALSE, type = 'lower',
         method = 'square', tl.cex = 0.8, tl.srt = 45, tl.col = 'red')





# UncertainORd

inverts <- read.csv(file = './Data/inverts_to_species.csv')
invert.uo <- ordinate_poisson(1000, inverts[,-1])



# bsmds

library(bsmds)
library(vegan)
inverts <- read.csv(file = './Data/inverts_to_species.csv')

source('./Code/bsmds_gp.R')

makedist.gp(thermometers2004, dist.fun = 'vegdist', dist.args = list(method = 'bray'))

m <-metaMDS(inverts[,-1], autotransform = FALSE)

nsim <- 9
out <- bsmds.gp(inverts[,-1], dist.fun = "vegdist", dist.data.arg = "x", dist.args=list(method="bray"), R=nsim, type="ordinal", iter.info=TRUE)

out <- bsmds.gp(thermometers2004, dist.fun = "vegdist", dist.data.arg = "x", dist.args=list(method="bray"), R=nsim, type="interval", iter.info=TRUE)

df <- do.call(rbind, out$X.i)
df <- data.frame(df)

df$grp <- rep(names(out$X.i), each = nsim)

library(tidyverse)
ggplot(df) + 
  geom_point(aes(x = X1, y = X2, col = grp), alpha = 0.5) +
  coord_equal()


#################  fitNMDS
library(fitNMDS)
# data(smoky)
x <- smoky$spe
b <- resamp_nmds(x, BS = nrow(x) * 0.8, B = 99, k=2) #BS is rows to keep, B = reps

df <- as.data.frame(b$points)
df.tf <- df %>% 
  rownames_to_column(var = 'tmp') %>%
  mutate(site = stringr::str_sub(tmp, 2)) %>%
  separate(site, into = c('site', 'rep')) %>%
  select(-c(tmp, rep))

# https://stackoverflow.com/questions/48690755/adding-convex-hull-to-ggplot-map  
hull.site <- df.tf %>%
  group_by(site) %>%
  slice(chull(V1, V2))


ggplot(df.tf) + 
  geom_point(aes(x = V1, y = V2, col = site), alpha = 0.5) +
  geom_polygon(data = hull.site, aes(V1, V2, col = site, fill = site), alpha = 0.05) +
  coord_equal()
                  



## gg stressplot
invert.dist <- vegdist(inverts)

invert.mds <- metaMDS(invert.dist, distance = 'bray', autotransform = FALSE, trace = 0)

object <- invert.mds
s <- stressplot(invert.mds)
g <- goodness(invert.mds)

s.df <- data.frame(s) %>%
  arrange(x)

plot(s$y ~ s$x)
lines(s.df$x, s.df$yf, type = "S", col = 'red')

ggplot(s.df) +
  geom_point(aes(x, y)) +
  geom_path(aes(x, yf), col = 'red', lwd = 1) +
  labs(x= 'Dissimilarity', y = 'Ordination distance') +
  theme_minimal()

rstress <- 1 - object$stress^2
ralscal <- cor(s$y, s$yf)^2

Rst <- format(rstress, digits = 3)
Ral <- format(ralscal, digits = 3)

stress <- sqrt(sum((s$y - s$yf)^2)/sum(s$y^2))
rstress <- 1 - stress





library(mvabund) ## Load a dataset from the mvabund package
data(spider)
y <- spider$abun
X <- scale(spider$x)
n <- nrow(y)
p <- ncol(y)

## NOTE: The values below MUST NOT be used in a real application;
## they are only used here to make the examples run quick!!!
example_mcmc_control <- list(n.burnin = 10, n.iteration = 100, 
                             n.thin = 1)

testpath <- file.path(tempdir(), "jagsboralmodel.txt")


## Example 1 - model with two latent variables, site effects, 
## 	and no environmental covariates
library(mvabund)
data(spider)
y <- spider$abun
X <- scale(spider$x)
n <- nrow(y)
p <- ncol(y)

library(boral)
testpath <- file.path(tempdir(), "jagsboralmodel.txt")
example_mcmc_control <- list(n.burnin = 10, n.iteration = 100, 
                             n.thin = 1)
spiderfit_nb <- boral(y, family = "negative.binomial", 
                      lv.control = list(num.lv = 2), row.eff = "fixed",  mcmc.control = example_mcmc_control, 
                      model.name = testpath)
spiderfit_nb$jags.model$BUGSoutput$sims.matrix

summary(spiderfit_nb)
mcmc.sims <- 