
mva.visualise <- function(X, which.plot = c('heat', 'single', 'warton'))
{
  require(ggplot2)
  
  X.long <- X %>%
    rownames_to_column(var = 'site') %>%
    pivot_longer(cols = -site, names_to = 'taxa', values_to = 'abund')
  
  
  if ('heat' %in% which.plot)
  {  
    heat.gg <- ggplot(X.long, aes(x = site, y = taxa)) + 
      geom_raster(aes(fill=log10(abund)), show.legend = FALSE) + 
      scale_fill_gradient(low="grey95", high="red") +
      labs(x = 'Site', y = 'Taxa') +
      theme_minimal() +
      theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
  }
##
  if ('single' %in% which.plot)
  {  
    n.sites <- data.frame(n = colSums(decostand(X, method = 'pa')))

    sites.hist.gg <- ggplot(n.sites) +
      geom_histogram(aes(x = n), binwidth = 1) +
      labs(x = 'No. of sites', y = 'Frequency') +
      theme_minimal()
  }

  ##
  if ('warton' %in% which.plot)
  {  
    warton <- data.frame(x = as.vector(as.matrix(X)),
                        y = rep(1:ncol(X), each = nrow(X)),
                        col = 'red')

    dot.warton.gg <- ggplot(data = warton) +
      geom_point(aes(x = log10(x+1), y  = y, col = col), alpha = 0.5, show.legend = FALSE) +
      labs(y= 'Taxa', x = 'Abundance (log10[x+1]') +
      theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
    theme_minimal()
  }
   
  plot.list <- vector(mode = 'list')
  heat.gg  
  plot.list[[1]] <- ifelse (exists('heat.gg'), heat.gg, NA)
  plot.list[[2]] <- ifelse (exists('sites.hist.gg'), sites.hist.gg, NA)
  plot.list[[3]] <- ifelse (exists('dot.warton.gg'), NA)

  plot.list  
}