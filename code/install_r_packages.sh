#!/bin/bash
# Install remaining necessary R packages:
Rscript -e "install.packages('https://cran.rstudio.com/src/contrib/calibrate_1.7.7.tar.gz', repos=NULL, type='source')"
Rscript -e "install.packages('https://cran.rstudio.com/src/contrib/msm_1.7.1.tar.gz', repos=NULL, type='source')"
Rscript -e "install.packages('https://cran.rstudio.com/src/contrib/RcppParallel_5.1.7.tar.gz', repos=NULL, type='source')"

# TODO: update links (e.g. https://cran.r-project.org/src/contrib/Archive/phytools/)
Rscript -e "install.packages('https://cran.rstudio.com/src/contrib/phytools_2.1-1.tar.gz', repos=NULL, type='source')"
