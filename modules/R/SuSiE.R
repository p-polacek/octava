install.packages("susieR")

library(susieR)
library(data.table)

setwd("C:/Users/polac001/OneDrive - Wageningen University & Research/PhD/HS/TBD_pipeline/modules/R")

ss <- fread("fm_test_tomato_MLM_num_locules.tsv", sep="\t")
#ss <- ss[ss$Chr == "SL4.0ch11", ]
#ss <- ss[ss$Pos >= 52752504, ]
#ss <- ss[ss$Pos <= 53752504, ]
ss <- ss[ss$Chr == "SL4.0ch02", ]
ss <- ss[ss$Pos >= 44189386, ]
ss <- ss[ss$Pos <= 46189386, ]

ss <- ss[ss$Pos != 44770232]
ss <- ss[ss$Pos != 45571700]
ss <- ss[ss$Pos != 45562846]

phe <- fread("../../data/Tomato/Agronomic.phe", sep="\t")$Number_of_Locules
#Z <- ss$Effect / ss$std.Error
r <- fread("num_locules_chr02_r.phased.vcor1", sep="\t")
N <- 335

#susie_plot(z_scores, y = "z")

fitted_rss1 <- susie_rss(#z=Z,
                         bhat=ss$Effect,
                         shat=ss$std.Error,
                         n=N,
                         R=r,
                         var_y=var(phe),
                         L=10,
                         estimate_residual_variance=FALSE,
                         max_iter=100)

summary(fitted_rss1)$cs
susie_plot(fitted_rss1, y="PIP", b=b)

ss$PIP <- fitted_rss1$pip
write.table(subset(ss, select=c("Chr", "Pos", "PIP")), file="num_locules_chr02_pip.txt", sep="\t", row.names=F)

### Example

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
