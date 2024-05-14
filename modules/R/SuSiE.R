#install.packages("optparse", repos="https://cloud.r-project.org")
#install.packages("susieR")

library(data.table)
library(optparse)
library(susieR)
library(glue)

#setwd("/mnt/c/Users/polac001/OneDrive - Wageningen University & Research/PhD/HS/TBD_pipeline/modules/R")
#setwd("C://Users/polac001/OneDrive - Wageningen University & Research/PhD/HS/TBD_pipeline/modules/R")

#setwd("/mnt/c/Users/polac001/OneDrive - Wageningen University & Research/PhD/HS/TBD_pipeline/")
#setwd("C://Users/polac001/OneDrive - Wageningen University & Research/PhD/HS/TBD_pipeline/")

option_list <- list(
  make_option(c("-a", "--assoc"),
              type="character",
              default=NULL,
              metavar="path_to_file",
              help=".tsv.gz file containing association analysis output data from GAPIT"),
  make_option(c("-r", "--r2"),
              type="character",
              default=NULL,
              metavar="path_to_file",
              help="tab-separated square matrix containng r2 for all markers in region."),
  make_option(c("-v", "--var"),
              type="character",
              default=NULL,
              metavar="variance",
              help="Variance of target phenotype values"),
  make_option(c("-n", "--num_phe"),
              type="character",
              default=NULL,
              metavar="num_phe",
              help="Number of phenotyped samples"),
  make_option(c("-c", "--chrom"),
              type="character",
              default=NULL,
              metavar="chromosome",
              help="Target chromosome."),
  make_option(c("-s", "--start"),
              type="integer",
              default=NULL,
              metavar="start",
              help="Target region at or past this position."),
  make_option(c("-e", "--end"),
              type="integer",
              default=NULL,
              metavar="end",
              help="Target region at or before this position."),
  make_option(c("-o", "--outdir"),
              type="character",
              default=NULL,
              metavar="path",
              help="Output directory.")
)

opt_parser = OptionParser(option_list=option_list,
                          epilogue="");
opts <- parse_args(opt_parser)
q()

if (0) {
  opts$assoc  <- "data/Demo/assoc/D/res/MLM.Plant_growth_habit.tsv.gz"
  opts$chrom <- "SL4.0ch06"
  opts$start <- 36133298
  opts$end   <- 47258699
  opts$var <- 4.349
  opts$num_phe <- 253
  opts$r2    <- "tmp/70aecc764e908d442f08e8202f3282ef1f5c81fe447f380c26d2996424c4c02e.phased.vcor2"
}

assoc <- fread(opts$assoc, sep="\t", check.names=F)
assoc <- assoc[assoc$Chr == opts$chrom, ]
assoc <- assoc[assoc$Pos >= opts$start, ]
assoc <- assoc[assoc$Pos <= opts$end, ]

#assoc <- assoc[assoc$Pos != 44770232]
#assoc <- assoc[assoc$Pos != 45571700]
#assoc <- assoc[assoc$Pos != 45562846]

r <- as.matrix(fread(opts$r2, sep="\t", check.names=F))
fitted_rss1 <- susie_rss(bhat=assoc$Effect,
                         shat=assoc[["std Error"]],
                         n=opts$num_phe,
                         R=r,
                         var_y=as.numeric(opts$var),
                         L=10,
                         estimate_residual_variance=FALSE,
                         max_iter=100)

#summary(fitted_rss1)$cs
#susie_plot(fitted_rss1, y="PIP", b=b)

assoc$PIP <- fitted_rss1$pip

setwd(opts$outdir)
write.table(subset(assoc, select=c("Chr", "Pos", "PIP")), file="pip.txt", sep="\t", row.names=F, quote=F)

### Example
q()

set.seed(1)

data(N3finemapping)
attach(N3finemapping)

n = nrow(X)
dim(Y)

#b <- true_coef[,1]
#plot(b, pch=16, ylab='effect size')
#which(b != 0)

sumstats <- univariate_regression(X, Y[,1])
z_scores <- sumstats$betahat / sumstats$sebetahat
susie_plot(z_scores, y = "z", b=b)

R <- cor(X)

fitted_rss1 <- susie_rss(bhat = sumstats$betahat,
                         shat = sumstats$sebetahat,
                         n = n,
                         R = R,
                         var_y = var(Y[,1]),
                         L = 10,
                         estimate_residual_variance = TRUE)
