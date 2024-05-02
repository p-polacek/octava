#install.packages("optparse", repos="https://cloud.r-project.org")
#install.packages("susieR")

library(data.table)
library(optparse)
library(susieR)
library(glue)

setwd("/mnt/c/Users/polac001/OneDrive - Wageningen University & Research/PhD/HS/TBD_pipeline/modules/R")
#setwd("C:/Users/polac001/OneDrive - Wageningen University & Research/PhD/HS/TBD_pipeline/modules/R")

setwd("/mnt/c/Users/polac001/OneDrive - Wageningen University & Research/PhD/HS/TBD_pipeline/")
#setwd("C:/Users/polac001/OneDrive - Wageningen University & Research/PhD/HS/TBD_pipeline/")

option_list <- list(
  make_option(c("-g", "--gwas"),
              type="character",
              default=NULL,
              metavar="path_to_file",
              help=".tsv.gz file containing gwas output data from GAPIT"),
  make_option(c("-r", "--r2"),
              type="character",
              default=NULL,
              metavar="path_to_file",
              help="tab-separated square matrix containng r2 for all markers in region."),
  make_option(c("-p", "--phenotype"),
              type="character",
              default=NULL,
              metavar="path_to_file",
              help=".csv file containing phenotype data, one accession per row and one trait per column, must have row and column names"),
  make_option(c("-t", "--trait"),
              type="character",
              default=NULL,
              metavar="trait",
              help="Trait to run analysis for"),
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

if (0) {
  opts$gwas  <- "data/Tomato/GWAS/All SNPs/assoc/MLM.Blossom_end_scar.tsv.gz"
  opts$chrom <- "SL4.0ch02"
  opts$start <- 41139896
  opts$end   <- 49323936
  opts$phenotype <- "data/Tomato/Agronomic.phe"
  opts$trait <- "Blossom_end_scar"
  opts$r2    <- "data/Tomato/GWAS/All SNPs/imp_filt/snp.ld.phased.vcor2"
}

gwas <- fread(opts$gwas, sep="\t")
gwas <- gwas[gwas$Chr == opts$chrom, ]
gwas <- gwas[gwas$Pos >= opts$start, ]
gwas <- gwas[gwas$Pos <= opts$end, ]

#gwas <- gwas[gwas$Pos != 44770232]
#gwas <- gwas[gwas$Pos != 45571700]
#gwas <- gwas[gwas$Pos != 45562846]

phe <- fread(opts$phenotype, sep="\t")[[opts$trait]]
phe <- phe[!is.na(phe)]

r <- as.matrix(fread(opts$r2, sep="\t"))
N <- length(phe)

fitted_rss1 <- susie_rss(bhat=gwas$Effect,
                         shat=gwas[["std Error"]],
                         n=N,
                         R=r,
                         var_y=var(phe),
                         L=10,
                         estimate_residual_variance=FALSE,
                         max_iter=100)

#summary(fitted_rss1)$cs
#susie_plot(fitted_rss1, y="PIP", b=b)

gwas$PIP <- fitted_rss1$pip

setwd(opts$outdir)
write.table(subset(gwas, select=c("Chr", "Pos", "PIP")), file="pip.txt", sep="\t", row.names=F, quote=F)

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
