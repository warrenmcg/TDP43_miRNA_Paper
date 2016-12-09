#!/usr/bin/Rscript

source("http://bioconductor.org/biocLite.R")
biocLite("ProMISe")
biocLite("DESeq2")
biocLite("SPIA")

install.packages("HiveR", repos="http://cran.us.r-project.org")
install.packages("gridExtra", repos="http://cran.us.r-project.org")
install.packages("optparse", repos="http://cran.us.r-project.org")