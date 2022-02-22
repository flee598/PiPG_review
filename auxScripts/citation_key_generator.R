# generate .bib file of packages, then extract the .bib key to be copied/pasted

library(bib2df)
library(tidyverse)

# bib file to create
fl <- "./package_references.bib"


# generate .bib for all packages used
knitr::write_bib(c("knitr", "rmarkdown", "pacman", "devtools", "tidyverse", "broom", "knitr", "janitor", "vegan", "analogue", "Rtsne", "umap", "ecotraj", "coenoflex", "indicspecies", "BiodiversityR", "ggfortify", "patchwork", "mdthemes", "ecole", "fitNMDS", "pairwiseAdonis", "ggordiplots", "ggvegan"), width = 60, file = fl)

# read in .bib
df <- bib2df(fl)

keys <- df %>%
  filter(CATEGORY == "MANUAL") %>%
  pull(BIBTEXKEY)


cat(paste0("@", keys), sep = ", ")
