mtx.descr <- function(X)
{
  X.pa <- vegan::decostand(X, method = 'pa')
  sing <- sum(colSums(X.pa) == 1) 
  pr.zeros <- 1 - (sum(X > 0) / prod(dim(X)))
  dbi <- mean(scale(X, center = FALSE, scale =  apply(X, 2, max)))
  
  list(n.elem = prod(dim(X)),
       n = dim(X)[1],
       p = dim(X)[2],
       singletons = round(sing, 3), 
       prop.zero = round(pr.zeros, 3), 
       dbi = round(dbi, 3))
}


plot.mds.gg <- function(mds, txt.col = 'blue', txt.x = 0, txt.y = 0, clusters = NULL, labels = FALSE, msl = 0.5, axes = c(1,2))
# msl: min seg length for ggrepel labels  
{
  require(tidyverse)
  
  xy <- data.frame(scores(mds)[,axes]) %>%
    rownames_to_column(var = 'site')
  n <- dim(scores(mds))[1]
  
  stress.lbl <- paste('Stress = ', round(mds$stress,3))
  
  if(is.null(clusters))
  {  
    plot.obj <- ggplot(data = xy) + 
      geom_point(aes(x = NMDS1, y = NMDS2), alpha  = 0.3) + 
      annotate(geom = 'text', x = txt.x, y = txt.y, label = stress.lbl, col = txt.col) +
      labs(x = 'NDMS 1', y = 'NMDS 2') +
      theme_minimal() +
      coord_equal() + 
      theme(legend.position = "none")
  } else 
  {    
    xy$cl <- clusters
    
    plot.obj <- ggplot(data = xy) + 
      geom_point(aes(x = NMDS1, y = NMDS2, col = cl), alpha = 0.5) + 
      annotate(geom = 'text', x = txt.x, y = txt.y, label = stress.lbl, col = txt.col) +
      labs(x = 'NDMS 1', y = 'NMDS 2') +
      theme_minimal() +
      coord_equal() 
  }
  
  if (labels == TRUE)
  {
    plot.obj <- plot.obj + 
      ggrepel::geom_text_repel(aes(x = NMDS1, y = NMDS2, label = site), min.segment.length = msl)
  }
    
  plot.obj
}  


plot.mds.bubble.gg <- function(mds, weights, rescale = 1, clusters = NULL, labels = FALSE,...)
# xy is the point location (e.g. the $points from MDS), weights the bubble size (e.g. from site-species mtx)
# rescale is the size of the bubbles (reduce if they're too big)
{
  xy <- data.frame(scores(mds)) %>%
    rownames_to_column(var = 'site') %>%
    mutate(cluster = clusters, weights = weights, plot_wgt = rescale * (weights / 4) + 1.5)
  
  n <- dim(scores(mds))[1]
  
  xy_pres <- filter(xy, weights > 0)
  xy_abs <- filter(xy, weights == 0)
  
  if(is.null(clusters))
  {  
    plot.obj <- ggplot() + 
      geom_point(data = xy_abs, aes(x = NMDS1, y = NMDS2), col = "grey", alpha = 0.4) + 
      geom_point(data = xy_pres, aes(x = NMDS1, y = NMDS2, size = plot_wgt), col = "black", pch = 1) + 
      # annotate(geom = 'text', x = txt.x, y = txt.y, label = stress.lbl, col = txt.col) +
      labs(x = 'NDMS 1', y = 'NMDS 2') +
      theme_minimal() +
      coord_equal() + 
      theme(legend.position = "none")
  } else 
  {    
    xy$cl <- as.factor(clusters)
    
    plot.obj <- ggplot(data = xy) + 
      geom_point(data = xy_abs, aes(x = NMDS1, y = NMDS2), alpha  = 0.4, col = "grey") + 
      
      geom_point(data = xy_pres, aes(x = NMDS1, y = NMDS2, size = plot_wgt, col = clusters), pch = 1) + 
      scale_color_distiller(type = "qual") +
      # annotate(geom = 'text', x = txt.x, y = txt.y, label = stress.lbl, col = txt.col) +
      labs(x = 'NDMS 1', y = 'NMDS 2') +
      theme_minimal() +
      coord_equal() 
  }
  
  if (labels == TRUE)
  {
    plot.obj <- plot.obj + 
      ggrepel::geom_text_repel(aes(x = NMDS1, y = NMDS2, label = site), min.segment.length = msl)
  }
  
  plot.obj
}

plot.envfit.gg <- function(mds, en, p.max = 0.05, txt.col = 'blue', txt.x = 0, txt.y = 0, clusters = NULL)
{
  require(ggplot2)
  
  xy <- data.frame(scores(mds))
  n <- dim(scores(mds))[1]
  
  stress.lbl <- paste('Stress = ', round(mds$stress,3))
  
  if(is.null(clusters))
  {  
    plot.obj <- ggplot(data = xy) + 
      geom_point(aes(x = NMDS1, y = NMDS2), alpha  = 0.8) + 
      # annotate(geom = 'text', x = txt.x, y = txt.y, label = stress.lbl, col = txt.col) +
      labs(x = 'NDMS 1', y = 'NMDS 2') +
      theme_minimal() +
      coord_equal() + 
      theme(legend.position = "none")
  } else 
  {    
    xy$cl <- as.factor(clusters)
    
    plot.obj <- ggplot(data = xy) + 
      geom_point(aes(x = NMDS1, y = NMDS2, col = cl), alpha  = 0.5) + 
      # annotate(geom = 'text', x = txt.x, y = txt.y, label = stress.lbl, col = txt.col) +
      labs(x = 'NDMS 1', y = 'NMDS 2') +
      theme_minimal() +
      coord_equal() 
  }
  
  # Build df for vectors and filter by p.values
  en.coord.cont <- data.frame(scores(en, "vectors") * ordiArrowMul(en)) 
  en.coord.cat <- data.frame(scores(en, "factors") * ordiArrowMul(en))

  # Add the vectors (code adapted from https://jkzorz.github.io/2020/04/04/NMDS-extras.html)
  if (nrow(en.coord.cont) > 0)
  {
    en.coord.cont <- en.coord.cont %>%
      rownames_to_column('variable') %>%
      mutate(p = en$vectors$pvals) %>%
      dplyr::filter(p <= p.max)

    plot.obj <- plot.obj +
      geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
                   data = en.coord.cont, alpha = 0.5, colour = "grey30", 
                   arrow = arrow(length = unit(0.2, "cm"))) +
      geom_text(data = en.coord.cont, aes(x = NMDS1, y = NMDS2, label = variable), colour = "grey30",
                fontface = "bold")
  }

  if (nrow(en.coord.cat) > 0)
  {
    en.coord.cat <- en.coord.cat %>%
      rownames_to_column('variable') %>%
      mutate(p = en$factors$pvals) %>%
      dplyr::filter(p <= p.max)

    plot.obj <- plot.obj +
      geom_point(data = en.coord.cat, aes(x = NMDS1, y = NMDS2),
                  shape = "diamond", size = 4, alpha = 0.6, colour = "navy") +
      geom_text(data = en.coord.cat, aes(x = NMDS1, y = NMDS2+0.04),
              label = row.names(en.coord.cat), colour = "navy", fontface = "bold")
  }
  
  plot.obj
}  

plot.shepp.gg <- function(dist, mds)
{
  s <- MASS::Shepard(dist, scores(mds), 2)
                          
  s.df <- data.frame(s) %>%
    arrange(x)
                          
  plot.obj <- ggplot(s.df) +
    geom_point(aes(x, y), alpha = 0.25) +
    geom_path(aes(x, yf), col = 'red', lwd = 1) +
    labs(x= 'Dissimilarity', y = 'Ordination distance') +
    theme_minimal()
  
  plot.obj
}  

plot.resamp <- function(b.df, legend = FALSE)
{
  require(ggplot2)
  
# Build convex hulls around the sites or clusters (sites by default unless clusters specified)
  hull.id <- b.df %>%
    group_by(site) %>%
    slice(chull(X1, X2))
  
# Plot
  mds.resamp.gg <- ggplot(b.df) + 
    geom_point(aes(x = X1, y = X2, col = site), alpha = 0.3, show.legend = FALSE) + 
    geom_polygon(data = hull.id, aes(X1, X2, col = site, fill = site, group = site), alpha = 0.05, show.legend = legend) +
    labs(x = 'NMDS1', y = 'NMDS2') +
    theme_minimal() +
    coord_equal()
  
  mds.resamp.gg
}


## plot a principle curve in PCA space using ggplot
# modified from https://github.com/gavinsimpson/analogue/blob/master/R/plot.prcurve.R
plot.prcurve.gg <- function(prc.obj, axes = 1:2, scaling = 0, segments = TRUE, site.label = TRUE)
{
  pred <- predict(prc.obj$ordination, prc.obj$s, type = "wa", scaling = scaling)[,axes]
  scrs <- scores(prc.obj$ordination, display = "sites", scaling = scaling, choices = axes)
  xlim <- range(scrs[,1], pred[,1])
  ylim <- range(scrs[,2], pred[,2])
  
  plot.obj <- ggplot() +
    geom_point(data = scrs, aes(x = PC1, y = PC2), size = 2, shape = 21) + 
    ggplot2:::limits(xlim, "x") +
    ggplot2:::limits(ylim, "y") +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_hline(yintercept = 0, linetype = 2) +
    coord_equal() +
    theme_minimal()
    
  
  if(segments) {
    segs <- data.frame(s = scrs, p = pred)
    plot.obj <- plot.obj + 
      geom_segment(data = segs, aes(x = s.PC1, y = s.PC2, xend = p.PC1, yend = p.PC2), col = "forestgreen") +
      geom_path(data = pred[prc.obj$tag, 1:2], aes(PC1, PC2), col = "red")
  }
  
  if(site.label) {
    scrs <- data.frame(scrs) %>%
      rownames_to_column(var = "site")
    
    plot.obj <- plot.obj + 
      ggrepel::geom_text_repel(data = scrs, aes(x = PC1, y = PC2, label = site))
    
  }
  
  plot.obj
}
# ggplot species response curves
# modified https://github.com/gavinsimpson/analogue/blob/master/R/plot.sppResponse.R

plot.sppResponse.gg <- function(spp.rsp, focal_spp = seq_along(spp.rsp), display = c("observed","fitted"))
{
  display <- match.arg(display)
  nams <- names(spp.rsp)

  ## process which - this could be logical
  if (is.logical(focal_spp))
    focal_spp <- which(focal_spp) ## yeah, really!
  
  plot.list <- list()
  idx <- 1
  
  for (i in focal_spp) {
    ox <- spp.rsp[[i]]$observed$gradient
    oy <- spp.rsp[[i]]$observed$response
    fx <- spp.rsp[[i]]$fitted.values$gradient
    fy <- spp.rsp[[i]]$fitted.values$response
    xlim <- range(ox, fx)
    ylim <- range(oy, fy)
    
    obs_dfr <- data.frame(ox, oy)
    fit_dfr <- data.frame(fx, fy)
      
    plot.list[[idx]] <- ggplot() +
        geom_point(data = obs_dfr, aes(ox, oy), shape = 21, size = 2) +
        geom_line(data = fit_dfr, aes(fx, fy), col = "red") +
        xlim(xlim) +
        ylim(ylim) +
        labs(x = "Gradient", y = "Abundance") +
        ggtitle(nams[i]) +
        theme_bw()
    
    idx <- idx + 1
    
  }
  plot.obj <- patchwork::wrap_plots(plot.list)
  plot.obj
}

## ggplot TSNE objects produced by Rtsne or UMAP via UMAP
plot.nonlin.gg <- function(obj, site_labels = NULL, labels = FALSE, clusters = NULL)
{
  
  if (class(obj) == "umap") { 
    if(is.null(site_labels)) {site_labels <- 1:nrow(obj$layout)}
    xy <- data.frame(site = site_labels, obj$layout)
    } else 
    {
      if(is.null(site_labels)) {site_labels <- 1:nrow(obj$Y)}
      xy <- data.frame(site = site_labels, obj$Y)
  }
  
  if (is_null(clusters)) {
    plot_obj <- ggplot(xy) +
      geom_point(aes(x = X1, y = X2 )) } else {
  
    plot_obj <- ggplot(xy) +
      geom_point(aes(x = X1, y = X2, col = factor(clusters)))    
      
    }
    
  if (labels & !is.null(site_labels)) {
    plot_obj <- plot_obj +
      ggrepel::geom_text_repel(aes(x = X1, y = X2, label = site))
    
  }  
  
  plot_obj + 
    labs( x = "Dimension 1", y = "Dimension 2") +
    coord_equal() +
    theme_bw() +
    theme(legend.position = "bottom")
}  


plot.coeno.gg <- function(x, f_grd = 1)
{
  
  i <- f_grd
  
  n_spp <- x$params$numspc
  grad_len <- x$params$grdlen[i]
  
  spp_curves <- matrix(nrow = grad_len + 1, ncol = n_spp)
  
  for (j in 1:n_spp) {
    val <- x$spcamp[j, i, ]
    spp_curves[,j] <- spp_amp_curve(val, grad_len, x$maxabu[j])
  }
  
  spp_curves_lg <- data.frame(spp_curves) %>%
    pivot_longer(cols = everything(), names_to = "species", values_to = "abund") %>%
    group_by(species) %>%
    mutate(grad_pos = 1:101) %>%
    arrange(species, grad_pos)
  
  plot.obj <- ggplot() +
    xlim(0, x$params$grdlen[i]) +
    ylim(0, max(x$maxabu)) +
    labs(x = "Gradient", y = "Abundance")
  
  
  plot.obj <- plot.obj +
    geom_line(data = spp_curves_lg, aes(x = grad_pos, y = abund, 
                                        col = species), show.legend = FALSE) +
    theme_bw()
  
  list(grad_plot = plot.obj, spp_curves = spp_curves_lg)
}

plot.coeno.rich.gg <- function(x)
{
  s <- rowSums(x$veg > 0)
  s_grad <- data.frame(grad_pos = x$site, s = s)
  rich_plot <- ggplot(s_grad, aes(x = grad_pos, y = s)) +
    geom_point(size = 3) +
    geom_smooth(method = "gam") +
    labs(x = "Gradient", y = "Richness") +
    theme_bw()
  
  list(s_grad = s_grad, rich_plot = rich_plot)
}

# Modified from the coenoflex::plot.coenoflex function
spp_amp_curve <- function(val, grad_len, max_abund)
{
  curve <- rep(0, grad_len + 1)
  
  for (k in 0:grad_len) {
    if (k < val[1])  curve[k + 1] <- 0
    else if (val[1] < k & k <= val[2]) 
      curve[k + 1] <- 2 * ((k - val[1])/(val[3] - val[1]))^2
    else if (val[2] < k & k <= val[3]) 
      curve[k + 1] <- 1 - 2 * ((val[3] - k)/(val[3] - val[1]))^2
    else if (val[3] < k & k <= val[4]) 
      curve[k + 1] <- 1 - 2 * ((k - val[3])/(val[5] - val[3]))^2
    else if (val[4] <- k & k <= val[5]) 
      curve[k + 1] <- 2 * ((val[5] - k)/(val[5] - val[3]))^2
    else curve[k + 1] <- 0
  }
  
  curve <- curve * max_abund
  curve        
}     