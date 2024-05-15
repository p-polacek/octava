#install.packages("optparse", repos="https://cloud.r-pro  ject.org")
#install.packages("BiocManager")
#BiocManager::install("snpStats")
#devtools::install_github("SFUStatgen/LDheatmap")
#devtools::install_github("jiabowang/GAPIT",force=TRUE)

library(jsonlite)
library(optparse)
library(GAPIT)
library(glue)

source("modules/R/save_optparse_params.R")
source("modules/R/Input_loader.R")

################################################################################
#### load and process input
################################################################################

if (!file.exists("modules/R/GAPIT_params.RDS")) {
  save_params()
}
option_list <- readRDS("modules/R/GAPIT_params.RDS")

opt_parser = OptionParser(option_list=option_list,
                          epilogue="");
opts <- parse_args(opt_parser)
input <- load_input(opts, "GAPIT")
models <- strsplit(opts$models, ",")[[1]]

################################################################################
#### Test input
################################################################################

if (0) {
  setwd("C://Users/polac001/OneDrive - Wageningen University & Research/PhD/HS/TBD_pipeline")
  opts$genotype <- "data/Demo/assoc/D/gt/gt.raw.gz"
  opts$phenotype <- "data/Demo/Agronomic.phe"
  opts$map <- "data/Demo/assoc/D/gt/gt.map.gz"
  opts$phenotype_tsv <- T
  opts$pop_structure <- "data/Demo/assoc/D/gt/gt.eigenvec"
  opts$q_tsv <- T
  opts$top_q <- 3
  opts$kinship <- "data/Demo/assoc/D/gt/gt.rel"
  opts$k_tsv <- T
  opts$outdir <- "data/Demo/assoc/D/res/"
  opts$model <- "MLM"
  opts$traits <- "Plant_growth_habit"
}

################################################################################
#### Save run parameters
################################################################################

setwd(opts$outdir)
for (model in models) {
  for (col in colnames(input$phe)[2:length(colnames(input$phe))]) {
    params <- c(
      (length(colnames(input$gt)) - 1),
      sum(!is.na(input$phe[[col]])),
      dim(input$Q)[2] - 1,
      "Yes", #is.null(input$K),
      var(input$phe[!is.na(input$phe[[col]]), col])
    )
    labels <- c(
      "Number of markers",
      "Number of samples",
      "Number of PCs",
      "Kinship",
      "Phenotype variance"
    )
    write(
      toJSON(list(labels, params), auto_unbox=TRUE),
      glue("{model}.{col}.params"))
  }
}
################################################################################
#### Run GWAS
################################################################################

dir.create("raw",
		   recursive = T,
		   showWarnings = F)
setwd("raw")

res <- GAPIT(Y=input$phe,
			 GD=input$gt,
			 GM=input$map,
			 KI=input$K,
			 CV=input$Q,
			 model=models
)
