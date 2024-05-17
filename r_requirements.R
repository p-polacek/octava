install.packages("optparse", repos="https://cloud.r-project.org")
install.packages("glue")
install.packages("jsonlite")
install.packages("data.table")
install.packages("R.utils")
install.packages("BiocManager")
install.packages("devtools")

library(BiocManager)
library(devtools)
BiocManager::install("snpStats")
devtools::install_github("SFUStatgen/LDheatmap")
devtools::install_github("jiabowang/GAPIT",force=TRUE)


install.packages("susieR")