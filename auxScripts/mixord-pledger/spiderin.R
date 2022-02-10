# spiderin.R

# -----------------------------------------------------------------
#             Define the data matrix
# -----------------------------------------------------------------

spider.df <- read.csv("spiderAbund.csv",T)

# 28 samples (rows), 12 species (columns).

spider.mat <- as.matrix(spider.df)

sitenames <- as.character(1:28)
spnames   <- colnames(spider.df)
# Shorten the names:
spnam <- c("A.acc","A.cun","A.fab","A.lut","A.per","A.alb",
           "P.lug","P.mon","P.nig","P.pul","T.ter","Z.spi")

dimnames(spider.mat) <- list(sitenames,spnam)
