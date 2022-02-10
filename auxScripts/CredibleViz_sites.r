CredibleViz_sites <- function (coord1, coord2, sites, type = "points", items = NULL)
{
  require(tidyverse)
  # hacked version of UncertainOrd::CredibleViz to olour by some site vector rather than by object
  # sites is a vector that contains the site codes
    num.pts <- ncol(coord1)
    num.sims <- nrow(coord1)
    
    coord1.thin <- coord1 %>% 
      as.data.frame() %>% 
      tidyr::gather() %>%
      dplyr::mutate(key = ordered(key, levels = unique(key)))
   
    coord2.thin <- coord2 %>% 
      as.data.frame() %>% 
      tidyr::gather() %>%
      dplyr::mutate(key = ordered(key, levels = unique(key)))
    
    combined <- data.frame(
        lv1 = coord1.thin$value,
        lv2 = coord2.thin$value,
        key = rep(1:num.pts, each = num.sims)
      )
    
    combined.summarized <- combined %>% 
      dplyr::group_by(key) %>%
      dplyr::summarise(mean.lv1 = mean(lv1), mean.lv2 = mean(lv2)) %>%
      dplyr::mutate(id = 1:num.pts) %>%
      bind_cols(sites = as.factor(sites))
    
    plot.obj <-
      ggplot2::ggplot(combined.summarized, ggplot2::aes(x = mean.lv1, y = mean.lv2)) + 
      #ggplot2::geom_text(ggplot2::aes(label = id, color = as.factor(id))) + 
      ggplot2::geom_text(ggplot2::aes(label = id, color = sites)) + 
      ggplot2::guides(color = FALSE) +
      ggplot2::xlab("latent variable one") + ggplot2::ylab("latent variable two") +
      ggplot2::theme_bw() + 
      ggplot2::ggtitle("") + 
      ggplot2::ylim(-3.1, 3.1) + 
      ggplot2::xlim(-3.1, 3.1)
    
    
    
    gg_color_hue <- function(n) {
      hues = seq(15, 375, length = n + 1)
      grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
    }
    
    n.sites <- length(unique(sites))
    cols = gg_color_hue(n.sites)
    
    
    if (type == "scatter") {
      for (item.numb in 1:length(items)) {
        sims.tmp <- combined %>% 
          dplyr::filter(key == items[item.numb]) %>%
          dplyr::mutate(sites = as.factor(sites[items[item.numb]]))
        
        plot.obj <- plot.obj + 
          ggplot2::geom_point(data = sims.tmp, ggplot2::aes(x = lv1, y = lv2, color = sites), alpha = 0.15)
        plot.obj
      }
    }
    
    else if (type == "circles") {
      for (item.numb in 1:length(items)) {
        sims.temp <- combined %>% 
          dplyr::filter(key == items[item.numb]) %>%
          dplyr::select(lv1, lv2) %>% 
          as.matrix
        
        fit <- MASS::cov.mve(sims.temp,
                        quantile.used = nrow(sims.temp) * 0.95,
                        nsamp = 100)
        
        points_in_ellipse <- sims.temp[fit$best,]
        ellipse_boundary <- stats::predict(cluster::ellipsoidhull(points_in_ellipse))
        ellipse_boundary_df <- data.frame(ellipse_boundary)
        ellipse_boundary_df$sites <- as.factor(sites[items[item.numb]])
        
        plot.obj <- plot.obj + 
          ggplot2::geom_path(
            data = ellipse_boundary_df,
            ggplot2::aes(x = V1, y = y, color = sites),
            linetype = 3)
      }
    }
    
    else if (type == "density") {
      for (item.numb in 1:length(items)) {
        sims.tmp <- combined %>% 
          dplyr::filter(key == items[item.numb])
        
        dens <- MASS::kde2d(sims.tmp[, 1], sims.tmp[, 2], n = 75)
        ord_d <- expand.grid(z1 = dens$x, z2 = dens$y) %>%
          tibble::as_tibble() %>% 
          dplyr::mutate(den = as.vector(dens$z),
                        sites = as.factor(sites[items[item.numb]]))
        
        plot.obj <- plot.obj + 
          ggplot2::geom_contour(data = ord_d, 
                                ggplot2::aes(x = z1, y = z2, z = den, color = sites))
        plot.obj
      }
    }
    return(list(plot.obj = plot.obj))
  }
