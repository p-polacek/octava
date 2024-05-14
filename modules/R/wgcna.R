#install.packages("WGCNA")
#install.packages("BiocManager")
#BiocManager::install("GO.db")
#BiocManager::install("impute")
#BiocManager::install("preprocessCore")

library(WGCNA)
library(data.table)

allowWGCNAThreads(4)

setwd("../TBD_pipeline/modules/R/")
data <- read.table(
  "../../data/Pepper v1.6 G2PSOL/all_quant.combined.quantiles.txt.gz",
  sep="\t",
  header=T,
  row.names=1)
data <- t(data[colnames(data)[4:length(colnames(data))]])

powers = c(c(1:10), seq(from = 12, to = 20, by = 2))
sft  <- pickSoftThreshold(data, powerVector = powers, verbose = 5) 
