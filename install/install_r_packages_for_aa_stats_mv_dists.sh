#!/bin/bash
# Install remaining necessary R packages.
Rscript -e "install.packages('https://cran.rstudio.com/src/contrib/calibrate_1.7.7.tar.gz', repos=NULL, type='source')"
Rscript -e "install.packages('https://cran.rstudio.com/src/contrib/Archive/msm/msm_1.7.1.tar.gz', repos=NULL, type='source')"
Rscript -e "install.packages('https://cran.rstudio.com/src/contrib/Archive/RcppParallel/RcppParallel_5.1.7.tar.gz', repos=NULL, type='source')"
Rscript -e "install.packages('https://cran.rstudio.com/src/contrib/Archive/phytools/phytools_2.1-1.tar.gz', repos=NULL, type='source')"
