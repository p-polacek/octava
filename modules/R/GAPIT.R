#install.packages("optparse", repos="https://cloud.r-project.org")
#install.packages("BiocManager")
#BiocManager::install("snpStats")
#devtools::install_github("SFUStatgen/LDheatmap")
#devtools::install_github("jiabowang/GAPIT",force=TRUE)

library(optparse)
library(glue)
library(GAPIT)

setwd("/mnt/c/Users/polac001/OneDrive - Wageningen University & Research/PhD/HS/TBD_pipeline/modules/R")
#setwd("C:/Users/polac001/OneDrive - Wageningen University & Research/PhD/HS/TBD_pipeline/modules/R")
#source("gapit_functions.R")
source("save_optparse_params.R")
source("Input_loader.R")

################################################################################
#### load and process input
################################################################################

if (!file.exists("GAPIT_params.RDS")) {
  save_params()
}
option_list <- readRDS("GAPIT_params.RDS")

setwd("/mnt/c/Users/polac001/OneDrive - Wageningen University & Research/PhD/HS/TBD_pipeline/")
#setwd("C:/Users/polac001/OneDrive - Wageningen University & Research/PhD/HS/TBD_pipeline/")
opt_parser = OptionParser(option_list=option_list,
                          epilogue="Created by Peter Polacek");
opts <- parse_args(opt_parser)

#opts <- list(genotype="data/Tomato/GWAS/snp-prune/imp_filt/imp.filt.raw",
#             phenotype="data/Tomato/Agronomic.phe",
#             map="data/Tomato/GWAS/snp-prune/imp_filt/imp.filt.map",
#             phenotype_tsv=T,
#             chrom="SL4.0ch02",
#             start=1000000,
#             end=2000000,
#             model="MLM",
#             traits="a")

input <- load_input(opts, "GAPIT")

models <- strsplit(opts$models, ",")[[1]]

################################################################################
#### Run GWAS
################################################################################

model_outdir <- file.path(opts$outdir, "raw")
dir.create(model_outdir,
		   recursive = T,
		   showWarnings = F)

setwd(model_outdir)

res <- GAPIT(Y=input$phe,
			 GD=input$gt,
			 GM=input$map,
			 KI=input$K,
			 CV=input$Q,
			 model=models
)

q()

######################## SNPs

setwd("/media/scratchpad_01/guest22/RP/Data/Acc_filtered/SNPs/")

gd <- fread("imp.het_filt.G2PSOLMerge.raw", sep="\t")
gd <- subset(gd, select=-c(FID, PAT, MAT, SEX, PHENOTYPE))

gm <- fread("imp.het_filt.G2PSOLMerge.map", header = F, sep="\t")
gm <- subset(gm, select=c(V2, V1, V4))


for (model in c("BLINK", "MLM")) {
  
  dir.create(file.path("/media/scratchpad_01/guest22/RP/Results/GAPIT_out/Metabolic/",
                       glue("SNPs-{model}")), showWarnings = F)
  setwd(glue("/media/scratchpad_01/guest22/RP/Results/GAPIT_out/Metabolic/SNPs-{model}/"))
  
  res <- GAPIT(Y=phe,
               GD=gd,
               GM=gm,
               KI=K,
               CV=Q,
               model=model
  )
}