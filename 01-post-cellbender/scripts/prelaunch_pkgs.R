# Pre-Launch package installation

# CRAN
cran.pkgs = c("tidyverse", "openxlsx", "compositions", "future", "cowplot", "ggbiplot", "devtools")
install.packages(cran.pkgs)

## We used Seurat version 4.2.1 for our analysis (started back in 2022). We have to manually select that version.
install.packages("https://cran.r-project.org/src/contrib/Archive/Seurat/Seurat_4.2.1.tar.gz", repos=NULL, type="source")

# Bioconductor
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

bioconductor.pkgs = c("glmGamPoi","DESeq2","destiny")
BiocManager::install(bioconductor.pkgs)

# Github
devtools::install_github('satijalab/seurat-wrappers')
devtools::install_github('msraredon/NICHES', ref = 'master')