knitr::opts_chunk$set(echo = TRUE)

# load pacman package + install if missing - used to load/install multiple packages at once
if (!require("pacman")) install.packages("pacman")

# load packages + install any that are missing
pacman::p_load(tidyverse, broom, knitr, janitor, ggfortify, patchwork,
               mdthemes, vegan, analogue, ggvegan, mvabund, devtools)

## # install packages from GitHub (not available on CRAN)
## devtools::install_github("phytomosaic/ecole")
## devtools::install_github("phytomosaic/fitNMDS")
## devtools::install_github("jfq3/ggordiplots")
## devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
## devtools::install_github("gavinsimpson/ggvegan")
## 

# load newly installed packages from GitHub
pacman::p_load(ecole, fitNMDS, ggordiplots, pairwiseAdonis, ggvegan) 

# load custom functions for plotting ordinations from vegan
source("auxScripts/helperPlots.r")

# Load data and tidy up
# invert abundances
inverts <- read.csv(file = "data/inverts_to_species.csv", row.names = 1)
# water chem at some sites
chem <- read.csv(file = "data/physicochemical.csv", row.names = 1)
# location data (WGS84 coords)
locations <- read.csv(file = "data/Site_Coords_WGS84.csv", row.names = 1)

invert_sites <- rownames(inverts)
chem_sites <- rownames(chem)

# standardise data
inverts_pa <- decostand(inverts, method = "pa")   # Presence-absence
inverts_l10 <- log10(inverts + 1)                 # log-10 transformed
inverts_hel <- decostand(inverts, method = "hel") # Hellinger transformed

# Get a subset of the water chemistry data 
chem <- janitor::clean_names(chem)

# Only keep sit4s with invert data and remove some less useful variables
chem_sub <- chem[chem_sites %in% invert_sites,]
chem_sub <- select(chem_sub, -order, -land_form, -taxa, -qmci,
                   -mci, -ept, -proportion_ept, -gambusia)

# Subset of chemical data for hypothesis testing
chem_sub_hyp <- chem_sub %>%
  select(do, temperature, velocity, proportion_silt, proportion_macrophyte,
         proportion_riparian_cover20, pasture, din)

# Binary vectors for riparian cover and native cover
rip_cover <- chem_sub$proportion_riparian_cover0_5 > 0.2
native_cover <- chem_sub$native > 0.2

# Load the Aotea veg data
# vegan uses row names, but the tidyverse doesn"t!
aotea <- read.csv(file = "Data/aotea_allyears_tidy.csv", row.names = 1)
aotea_pa <- decostand(aotea, method = "pa")           # convert to presence-absence
aotea_pa_nosing <- aotea_pa[,colSums(aotea_pa) > 1]   # PA without any singletons

# Make the site labels more useful
dfr <- data.frame(site = rownames(aotea_pa))
aotea_site <- separate(dfr, site, into = "code", "_") # ignore the warning
aotea_site$code <- str_sub(aotea_site$code, end = 2)
rm(dfr)

# First, lengthen the data to make it easier for ggplot
aotea.long <- aotea %>%
  rownames_to_column(var = "site") %>%
  pivot_longer(cols = -site, names_to = "taxa", values_to = "abund")

inverts.long <- inverts %>% 
  rownames_to_column(var = "site") %>%
  pivot_longer(cols = -site, names_to = "taxa", values_to = "abund")

# Histogram of number of sites each spp occurs at for Lee et al.
n.sites.iv <- data.frame(n = colSums(inverts > 0))

sites.hist.iv.gg <- ggplot(n.sites.iv) + 
  geom_histogram(aes(x = n), binwidth = 1) +
  labs(x = "No. of sites", y = "Frequency") +
  theme_minimal()

# A heat map of abundance for species x site for Lee et al. data
inverts.long$taxa_num <- as.numeric(factor(inverts.long$taxa,
                                           levels = unique(inverts.long$taxa)))

heat.iv.gg <- ggplot(inverts.long, aes(x = site, y = taxa_num)) +
  geom_raster(aes(fill=log10(abund)), show.legend = FALSE) +
  scale_fill_gradient(low="grey99", high="red") +
  labs(x = "Site", y = "Taxa") +
  theme_minimal() + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# A dotplot to replicate that in Wang et al. 2012 Methods Ecol Evol
warton.iv <- data.frame(x = as.vector(as.matrix(inverts)), 
                        y = rep(1:77, each = 29), 
                        col = chem_sub$native > 0.4)

dot.wang.iv.gg <- ggplot(data = warton.iv) + 
  geom_point(aes(x = log10(x+1), y = y, col = col),
             alpha = 0.8, size = 2, show.legend = FALSE) + 
  scale_color_brewer(type = "qual") +
  labs(y= "Taxa", x = "Abundance (log<sub>10</sub>[x+1])") +
  theme(axis.text.x=element_blank(), 
          axis.ticks.x=element_blank()) +
  md_theme_minimal()   # mdthemes lets us use markdown in graph notations

# Bind them together with patchdown
fig1 <- sites.hist.iv.gg + heat.iv.gg + dot.wang.iv.gg +
plot_annotation(tag_levels = "a")

fig1

# freshwater invert data
# calculate mean and variance
inverts_mv <- data.frame(mean = colMeans(inverts), 
                         var = apply(inverts , 2, var))

# plot
invert.mean.gg <- ggplot(data = inverts_mv) + 
  geom_point(aes(x = log10(mean), y = log10(var)), alpha = 0.8, size = 2) + 
  labs(y = "Variance (log<sub>10</sub> scale)",
       x = "Mean (log<sub>10</sub> scale)") +
  md_theme_minimal() 

# Aotea vegetation data
# calculate mean and varaince
aotea_mv <- data.frame(mean = colMeans(aotea_pa), 
                       var = apply(aotea_pa , 2, var))

# plot
veg.mean.gg <- ggplot(data = aotea_mv) + 
  geom_point(aes(x = log10(mean), y = log10(var)), alpha = 0.8, size = 2) + 
  labs(y = "Variance (log<sub>10</sub> scale)",
       x = "Mean (log<sub>10</sub> scale)") +
  md_theme_minimal() 

# combine plots
fig2 <- invert.mean.gg + veg.mean.gg +
  plot_annotation(tag_levels = "a")

fig2

# Non-weighted using the Bray_Curtis
pcoa_iv <- wcmdscale(d = vegdist(inverts), eig = TRUE, add = "lingoes")

# get the weights (row sums)
wgt_iv <- iv <- rowSums(inverts)/sum(inverts)

 # Weighted
pcoa_wgt_iv <- wcmdscale(d = vegdist(inverts), w = wgt_iv, eig = TRUE, add = "lingoes")

# Get the coordinates (using vegan::scores)
pcoa_scores <- data.frame(scores(pcoa_iv)[,1:2],
                          scores(pcoa_wgt_iv)[,1:2],
                          site = rownames(inverts),
                          w = wgt_iv) %>%
  rename(Dim1_wgt = Dim1.1, Dim2_wgt = Dim2.1)

# Get the eigenvalues and their proportional value
pcoa_eig <- data.frame(eig = pcoa_iv$eig[1:10], 
                       expl = round(pcoa_iv$eig[1:10] / sum(pcoa_iv$eig[1:10]) * 100, 2))

pcoa_wgt_eig <- data.frame(eig = pcoa_wgt_iv$eig[1:10], 
                           expl = round(pcoa_wgt_iv$eig[1:10] / sum(pcoa_wgt_iv$eig[1:10]) * 100, 2))

# Make some axis labels
xl <- paste("Dimension 1 (",pcoa_eig$expl[1], "%)", sep = "")
yl <- paste("Dimension 2 (", pcoa_eig$expl[2], "%)", sep = "")

xl_w <- paste("Dimension 1 (",pcoa_wgt_eig$expl[1], "%)", sep = "")
yl_w <- paste("Dimension 2 (", pcoa_wgt_eig$expl[2], "%)", sep = "")

## unweighted plot
pcoa.gg <- ggplot(pcoa_scores) + 
  geom_text(aes(x = Dim1, y = Dim2, label = site), size = 3) + 
  labs(x = xl,
       y = yl) + 
  coord_equal() + 
  theme_bw()

# eigenvalues plot
pcoa.eig.gg <- ggplot(pcoa_eig) + 
  geom_bar(aes(x = as.factor(1:10), y = expl),
           stat = "identity") + 
  labs(x = "Dimension",
       y = "Eigenvalue (% total)") +
  theme_bw()

## weighted plot
pcoa.wgt.gg <- ggplot(pcoa_scores) + 
  geom_text(aes(x = Dim1_wgt, y = Dim2_wgt, label = site, size = w)) + 
  scale_size_continuous(range  = c(1, 4)) +
  labs(x = xl_w,
       y = yl_w,
       size = "weight") + 
  coord_equal() + 
  theme_bw()

# eigenvalues plot
pcoa.eig.wgt.gg <- ggplot(pcoa_wgt_eig) + 
  geom_bar(aes(x = as.factor(1:10), y = expl), stat = "identity") + 
  labs(x = "Dimension",
       y = "Eigenvalue (% total)") +
  theme_bw()

# gather plots
fig3 <- pcoa.gg + pcoa.eig.gg + pcoa.wgt.gg + pcoa.eig.wgt.gg + 
    plot_layout(widths = c(1, 1)) + plot_annotation(tag_levels = 'a') & theme(legend.position = 'bottom')


fig3

# nMDS of the Lee et al. invertebrate data

# Bray-curtis distance matrix
invert_dist <- vegdist(inverts)

# Turn off automatic rescaling. Trace controls output verbosity
invert_mds <- metaMDS(invert_dist, autotransform = FALSE, trace = 1)

# Make the plots (using custom plot.mds.gg function from helperPlots.r)
iv.mds.gg <- plot.mds.gg(invert_mds, txt.x = -0.35,
                         txt.y = -0.35,labels = TRUE,
                         msl = 0.25)


# NOT USED - remove? 

#nMDS coloured by native veg cover
# iv.mds_cluster.gg <- plot.mds.gg(invert_mds,
#                                  txt.x = -0.4, txt.y = -0.4,
#                                  clusters = native_cover) 

# Shepperd plots (distance in nMDS vs distance in dissimilarity matrix)
# iv.shep.gg <- plot.shepp.gg(invert_dist, invert_mds)



# nMDS with points scaled in size by abundance of _Physa_
wgt <- inverts$Physa / 10
iv.bubble.gg <- plot.mds.bubble.gg(invert_mds, weights = wgt)

# Fit the chemical data as environmental vectors (envfit) and plot
iv_chem_fit <- envfit(invert_mds, chem_sub_hyp)
iv.fit.gg <- plot.envfit.gg(invert_mds, iv_chem_fit, p.max = .2,
                            clusters = native_cover)

# Or we could synthesis all the chem data with a PCA and fit those axes
chem_pca <- princomp(chem_sub, cor = TRUE) 
chem_pca_scores <- scores(chem_pca)
iv_pca_fit <- envfit(invert_mds, chem_pca_scores[,1:8])

iv.fitpca.gg <- plot.envfit.gg(invert_mds,
                               iv_pca_fit,
                               p.max = .1,
                               clusters = native_cover)

# gather plots
fig4 <- (iv.mds.gg + iv.bubble.gg + iv.fit.gg + iv.fitpca.gg) +
  plot_annotation(tag_levels = "a") + 
  plot_layout(ncol = 2, guides ="collect") & 
  theme(legend.position = "bottom")

fig4

par(mfrow = c(2,2))

din_surf_k5 <- ordisurf(x = invert_mds, y = chem_sub$din,
                        knots = 5, main = "DIN, k = 5")

pca_surf_k5 <- ordisurf(x = invert_mds, y = chem_pca_scores[,1],
                        knots = 5, main = "PCA comp 1, k = 5")


din_surf_k10 <- ordisurf(x = invert_mds, y = chem_sub$din,
                         knots = 10, main = "DIN, k = 10")

pca_surf_k10 <- ordisurf(x = invert_mds, y = chem_pca_scores[,1],
                         knots = 10, main = "PCA comp 1, k = 10")

# nMDS of the Aotea veg data (presence-absence)
aotea_dist <- vegdist(aotea_pa)

# This code looks for a saved version of the nMDS and uses it if available (assumed stored in a folder called rds) - this speeds the analysis, which can be slow.
if (file.exists("rds/aotea_mds.rds")) { 
  aotea_mds <- readRDS("rds/aotea_mds.rds") } else {
  n_cores <- parallel::detectCores() - 1  
  aotea_mds <- metaMDS(aotea_dist, autotransform = FALSE, trace = 1,
                       parallel = n_cores)
  aotea_mds_3d <- metaMDS(aotea_dist, autotransform = FALSE,
                          try = 100, trace = 1, k = 3,
                          parallel = n_cores)
  
  saveRDS(aotea_mds, "rds/aotea_mds.rds") 
}

# plot Shepperd plot
aotea.shep.gg <- plot.shepp.gg(aotea_dist, aotea_mds)

# plot ordination
aotea.mds.gg <- plot.mds.gg(aotea_mds, txt.x = -0.6, txt.y = -0.35,
                            clusters = factor(aotea_site$code)) +
                theme(legend.position = "bottom")

# nMDS with points scaled to the abundance of the late successional tree, Beilschmiedia tawa
aotea.bubble <- plot.mds.bubble.gg(aotea_mds, weights = aotea$Beitaw)

# gather plots
fig5 <- (aotea.shep.gg + aotea.mds.gg) + 
  plot_annotation(tag_levels = "a")



fig5

aotea.bubble


# Fit PRC here with the complexity fixed across all species (vary)
# The method here specifies the starting conditions (first axis of a method)
# Here we use a GAM as the smoother (slower than spline)
inverts_pc <- prcurve(inverts_l10, method = "ca", trace = TRUE,
                      vary = TRUE, smoother = smoothGAM,
                      penalty = 1.4)
inverts_pc

plot.prcurve.gg(inverts_pc) +
  ggtitle("Final solution for principle curve")

# taxa responses
invert_responses <- sppResponse(inverts_pc)

# subset common taxa (> 60% sites)
focal_spp <- chooseTaxa(inverts_l10, n.occ = 17, value = FALSE)

# plot
plot.sppResponse.gg(invert_responses, focal_spp = focal_spp)


# load data
data(tikus, package = 'mvabund')

# log transform
tikus_log <- log(tikus$abund + 1)
tikus_survey <- tikus$x$time

# Unconstrained ordination of the Tikus coral reef data
# nMDS - increased try and maxit for convergence
tikus_nmds <- metaMDS(tikus_log, autotransform = TRUE,
                      wascores = FALSE,
                      try = 100, maxit = 500, trace = 1)

# plot
tik.mds.gg <- ggplot(data = as.data.frame(scores(tikus_nmds))) +
  geom_point(aes(NMDS1, NMDS2, col = as.factor(tikus$x$time)), size = 3) +
  scale_color_brewer(name = "Year", palette = "PuRd") +
  coord_equal() +
  theme_minimal()

# Constrained (by time) ordination of the Tikus coral reef data
# Note use of formula notation
tikus_rd <- dbrda(tikus_log ~ tikus_survey, distance = "bray")
tikus_rd_scores <- data.frame(scores(tikus_rd)$sites)

#plot
tik.rd.gg <- ggplot(data = tikus_rd_scores) + 
  geom_point(aes(dbRDA1, dbRDA2, col = as.factor(tikus_survey)), size = 3) +
  scale_color_brewer(name = "Year", palette = "PuRd") +
  coord_equal() +
  theme_minimal()

fig6 <- tik.mds.gg + tik.rd.gg + 
    plot_annotation(tag_levels = "a") +
    plot_layout(guides = "collect", widths = c(1, 1)) &
    theme(legend.position = "bottom")

fig6

invert_rda <- rda(inverts_hel ~ ., data = chem_sub_hyp) # RDA on subset of env data
invert_rda_r2 <- round(RsquareAdj(invert_rda)$adj.r.squared, 3) # adjusted r2 (via the vegan package)

# Plot using ggvegan (as an example)  
rda.iv.gg <- autoplot(invert_rda, layers = c("biplot", "sites"), geom = "text") +
  theme_minimal() + 
  theme(legend.position = "none")

rda.iv.gg

invert_cca_all <- cca(inverts_l10, chem_sub)
invert_cca_red <- cca(inverts_l10, chem_sub_hyp)


# Plot them - autoplot here is from ggvegan
invert.cca.all.gg <- autoplot(invert_cca_all, layers = c("biplot","sites"), geom = "text") + 
  theme_bw() + 
  theme(legend.position = "none")

invert.cca.red.gg <- autoplot(invert_cca_red, layers = c("biplot","sites"), geom = "text") + 
  theme_bw() + 
  theme_bw() + theme(legend.position = "none")

(invert.cca.all.gg + invert.cca.red.gg) + 
  plot_annotation(tag_levels = "a")


invert_dbrda <- capscale(inverts ~ ., chem_sub_hyp, dist="bray")
dbrda.invert.gg <- autoplot(invert_dbrda, layers = c('sites', 'biplot')) +
  theme_bw() + 
  theme(legend.position = 'none')

dbrda.invert.gg

full_mod_dbrda <- capscale(inverts_hel ~ ., data = chem_sub_hyp, dist = "bray")  # full model
null_mod_dbrda <- capscale(inverts_hel ~ 1, data = chem_sub_hyp, dist = "bray")  # intercept-only model

sel_mod_dbrda <- ordistep(null_mod_dbrda, 
                          scope = formula(full_mod_dbrda), 
                          direction = "both", 
                          permutations = how(nperm = 499), trace = FALSE ) 


anova(sel_mod_dbrda)

full_mod_dbrda_r2 <- round(RsquareAdj(full_mod_dbrda)$adj.r.squared, 3)
sel_mod_dbrda_r2 <- round(RsquareAdj(sel_mod_dbrda)$adj.r.squared, 3)


dbrda.full.gg <- autoplot(full_mod_dbrda, layers = c("sites", "biplot"), geom = "text") +
  ggtitle(paste("Full model. Adjusted R2 = ", full_mod_dbrda_r2)) +
  theme_minimal()

dbrda.sel.gg <- autoplot(sel_mod_dbrda, layers = c("sites", "biplot"), geom = "text") +
  ggtitle(paste("Reduced model. Adjusted R2 = ", sel_mod_dbrda_r2)) +
  theme_minimal()

dbrda.full.gg + dbrda.sel.gg + 
  plot_layout(guides = "collect")

full_mod <- rda(inverts_hel ~ ., chem_sub_hyp) # Model with all explanatory variables

mod_aov <- anova(full_mod)  # This is a test for the overall model
term_aov <- anova(full_mod, by = "term") # can also do marginal tests - see vegan help
axis_aov <- anova(full_mod, by = "axis")

mod_aov  # Model is significant
term_aov # look for significant terms  
axis_aov # and axes



# Define number of simulations to do - slow for large data but could parallelise!
# We'll use fitNDMS::resamp_ndms to do this.  BS is the bootstrap size and B the no., of sims
# Underneath it is using vegan::metaMDS

nsims_iv <- 199 
nsims_aotea <- 5  # need more for analysis but this is **slow**! 

# Because this is slow we'll cache the op in a directory called 'rds'
  if (file.exists("rds/invert_resamp.rds")) {
      b_iv <- readRDS("rds/invert_resamp.rds") } else { 
      set.seed(1767283)    # set seed for reproducibility
      b_iv <- resamp_nmds(inverts, BS = nrow(inverts) * 0.8, B = nsims_iv, k = 2)
      saveRDS(b_iv, "rds/invert_resamp.rds") 
  }
  
  # as per above need to check convergence...
  if (file.exists("rds/aotea_resamp.rds")) { 
      b_aotea <- readRDS("rds/aotea_resamp.rds") } else { 
      set.seed(4337659)
      b_aotea <- resamp_nmds(aotea_pa, BS = nrow(aotea_pa) * 0.8, B = nsims_aotea, k = 2)
      saveRDS(b_aotea, "rds/aotea_resamp.rds") 
      }

# Code here will trigger some warnings re expected pieces and coercion - can ignore
# Note that the number of times a site is sampled will not necessarily equal n_sims 
# It is c. 0.8 x sims

  b_iv_dfr <- data.frame(b_iv$points) %>% 
    rownames_to_column(var = 'site') %>% 
    arrange(site) %>% 
    separate(site, into = c('site', 'rep'), sep = "\\.") %>% 
    mutate(rep = as.numeric(rep)) %>% 
    mutate(rep = if_else(is.na(rep), 1 , as.numeric(rep) + 1))

  # TO DO - this might need checking
  b_aotea_dfr <- data.frame(b_aotea$points) %>% 
    rownames_to_column(var = 'pcq.site') %>% 
    arrange(pcq.site) %>% 
    separate(pcq.site, into = c('pcq_point', 'rep'), sep = "\\.") %>%
    mutate(site = gsub("_.*","", pcq_point)) %>%
    mutate(rep = if_else(is.na(rep), 1 , as.numeric(rep) + 1))
  
resamp.iv.gg <- plot.resamp(b_iv_dfr) 
resamp.aotea.gg <- plot.resamp(b_aotea_dfr, legend = TRUE)

resamp.iv.stress.gg <- ggplot(data = data.frame(s = b_iv$stresses)) +
  geom_histogram(aes(s), bins = 20) + labs(x = 'Stress', y = 'Frequency') +
  theme_minimal()

resamp.aotea.stress.gg <- ggplot(data = data.frame(s = b_aotea$stresses)) + 
  geom_histogram(aes(s), bins = 10) + 
  labs(x = 'Stress', y = 'Frequency') + 
  theme_minimal()

(resamp.iv.gg + resamp.iv.stress.gg) / (resamp.aotea.gg + resamp.aotea.stress.gg) + 
  plot_annotation(tag_level = 'a') +
  plot_layout(guides = "collect") & 
  theme(legend.position = 'bottom')


n_perm <- 499     # reduce if you want to speed up!
n_cores <- parallel::detectCores() - 1      # cores for parallelisation (set to 1 for none)

aotea_dist <- vegdist(aotea_pa, method = "bray")

# ANOSIM on 'site
aotea_anosim <-anosim(aotea_dist, aotea_site$code, permutations = n_perm, parallel = n_cores)  

# PERMANOVA on 'site'
aotea_adonis <- adonis2(aotea_dist ~ aotea_site$code, permutations = n_perm, parallel = n_cores)  

# You can't do pairwise comparisons in adonis2 but the pairwiseAdonis pakcage allows this
### TO DO

# Test for homogeneity of variance (scatter)
aotea_mod <- betadisper(aotea_dist, aotea_site$code)
aotea_mod_perm <- permutest(aotea_mod, pairwise = TRUE, permutations = n_perm)
aotea_mod_HSD <- TukeyHSD(aotea_mod)

# permuted stats
aotea_perm_dfr <- data.frame(perm = c(permustats(aotea_anosim)$permutations, 
                                permustats(aotea_adonis)$permutations),
                                test = rep(c('ANOSIM', 'PERMANOVA'), each = 499))

# observed stat (for comparison)
aotea_stat_dfr <- data.frame(stat = c(permustats(aotea_anosim)$statistic, 
                                      permustats(aotea_adonis)$statistc),
                            test = c('ANOSIM', 'PERMANOVA'))

aotea.sim.gg <- ggplot(aotea_perm_dfr) + 
  geom_histogram(aes(x = perm)) +
  geom_vline(data = aotea_stat_dfr, aes(xintercept = stat), col = 'red') +
  facet_wrap(~test, scales = 'free_x') + 
  labs (x = 'Permutation statistics', y = 'Frequency') + 
  theme_minimal()


aotea.sim.gg


aotea.adonis.td <- tidy(aotea_adonis)
disp.rc.plot <- plot(aotea_mod, main = "")

inverts_dispw <- dispweight(inverts)   # vegan::dispweight implements the method

iv_var <- apply(inverts, 2, var) 
iv_mean <- apply(inverts, 2, mean)
iv_sig <- as.factor(ifelse(attr(inverts_dispw, "p") < 0.01, 1, 2))
iv_dfr <- data.frame(m = iv_mean, v = iv_var, s = iv_sig)

mv.raw.iv <- ggplot(iv_dfr) + 
  geom_point(aes(x = m , y = v, size = m /v, col = s)) + 
  geom_abline(slope = 1, intercept = 0) + 
  scale_x_log10() + 
  scale_y_log10() + 
  labs(x = "Mean abundance", y = "Var. abundance") +
  theme_bw()

invert_wgt_mds <- metaMDS(inverts_dispw, distance = 'bray', autotransform = FALSE, trace = 1)

iv.mds_wgt.gg <- plot.mds.gg(invert_wgt_mds, txt.x = -0.9, txt.y = -0.8, labels = TRUE)

(mv.raw.iv)
(iv.mds.gg + iv.mds_wgt.gg)

simper_veg <- simper(aotea_pa, group = aotea_site$code)
# simper_veg


if (file.exists("rds/aotea_simper_perm.rds")) { 
  aotea_mds <- readRDS("rds/aotea_simper_perm.rds") } else {
  n_cores <- parallel::detectCores() - 1  
  n_perm <-19
  simper_perm_veg <- simper(aotea_pa, group = aotea_site$code, permutations = n_perm, parallel = n_cores, trace = FALSE)
  
  saveRDS(aotea_mds, "rds/aotea_simper_perm.rds") 
  }


# summary(simper_perm_veg)

library(indicspecies)
n_perm <- 99
aotea_indic <- multipatt(aotea_pa, aotea_site$code, control = how(nperm = n_perm))
summary(aotea_indic, indvalcom = TRUE, alpha = 0.01) # summarise, showing A and B for species with p <= .01


awana_indicators <- indicators(aotea_pa, cluster = aotea_site$code, group = "AW", At = 0.7, Bt = 0.3, max.order = 2, verbose = TRUE) # set verbose FALSE to turn off the output

summary(awana_indicators)
print(awana_indicators)

n_boot <- 19
awana_indicators_boot <- indicators(aotea_pa, cluster = aotea_site$code, group = "AW", At = 0.7, Bt = 0.3, max.order = 2, nboot.ci = n_boot, verbose = FALSE) # set verbose FALSE to turn off the output

summary(awana_indicators_boot)
print(awana_indicators_boot$A) # check out the A values
