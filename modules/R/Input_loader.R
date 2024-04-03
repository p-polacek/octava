#install.packages("optparse", repos="https://cloud.r-project.org")

library(optparse)
library(data.table)
library(dplyr)

################################################################################
#### Validate input
################################################################################

load_input <- function(opts, tool) {

    if (is.null(opts$phenotype)) {
      stop("Phenotype file must be provided.", call.=FALSE)
    }

    if (is.null(opts$genotype)) {
      stop("Genotype file must be provided.", call.=FALSE)
    }

    if (is.null(opts$map)) {
      stop("Map file must be provided.", call.=FALSE)
    }

    if (!file.exists(opts$phenotype)) {
      stop(glue("Phenotype file {opts$phenotype} not found."), call.=FALSE)
    }

    if (!file.exists(opts$genotype)) {
      stop(glue("Genotype file {opts$genotype} not found."), call.=FALSE)
    }

    if (!file.exists(opts$map)) {
      stop(glue("Map file {opts$map} not found."), call.=FALSE)
    }

    if (!is.null(opts$pop_structure)) {
      if (!file.exists(opts$pop_structure)) {
        stop(glue("File {opts$pop_structure} not found."), call.=FALSE)
      }
      Q_given <- TRUE
    } else {
      Q_given <- FALSE
      Q <- NULL
    }

    if (!is.null(opts$kinship)) {
      if (!file.exists(opts$kinship)) {
        stop(glue("File {opts$kinship} not found."), call.=FALSE)
      }
      K_given <- TRUE
    } else {
      K_given <- FALSE
      K <- NULL
    }

################################################################################
#### Load input files
################################################################################

    phe_sep <- if (opts$phenotype_tsv) "\t"  else ","
    phe <- read.csv(opts$phenotype,
                    sep=phe_sep,
                    header=T)

    if (length(colnames(phe)) == 0) {
      stop("No phenotypes found or wrong separator specified.",
           call.=FALSE)
    }

    gt <- fread(opts$genotype,
                sep="\t")

    if (length(colnames(gt)) == 1) {
      stop("Wrong separator specified.",
           call.=FALSE)
    }
    if (length(colnames(gt)) == 6) {
      stop("No genotypes found.",
           call.=FALSE)
    }
	gt <- subset(gt, select=-c(FID, PAT, MAT, SEX, PHENOTYPE))

	map <- fread(opts$map,
  #map <- read.csv(opts$map,
  #map <- read.table(opts$map,
                    sep="\t",
                    header=F)

    if (length(colnames(map)) != 4) {
      stop("Map file does not contain all columns or wrong separator specified.",
           call.=FALSE)
    }
	#colnames(map) <- c("CHR", "ID", "DIST", "POS")
	map <- subset(map, select=c(V2, V1, V4))

    if (length(rownames(map)) != length(colnames(gt)) - 1) {
      stop("Genotype and map files have different number of markers.",
           call.=FALSE)
    }

    if (Q_given) {
      q_sep <- if (opts$q_tsv) "\t"  else ","
      Q <- read.csv(opts$pop_structure,
                    sep=q_sep,
                    header=T,
                    row.names=1)

      if (length(colnames(Q)) == 0) {
        stop("No data found or wrong separator specified for Q file.",
             call.=FALSE)
      }
    }

    if (K_given) {
      k_sep <- if (opts$k_tsv) "\t"  else ","
      K <- read.csv(opts$kinship,
                    sep=k_sep,
                    header=T,
                    row.names=1)

      if (length(colnames(K)) == 1) {
        stop("No data found or wrong separator specified for K file.",
             call.=FALSE)
      }
    }

################################################################################
#### Filter input files
################################################################################

	phe_accessions <- t(phe[, 1])
	gt_accessions  <- t(gt[, 1])
  filtered_accessions <- intersect(phe_accessions, gt_accessions)

  phe <- phe[phe_accessions %in% filtered_accessions,]
	gt <- gt[gt_accessions %in% filtered_accessions,]

  columns <- c(colnames(phe)[1], strsplit(opts$traits, ",")[[1]])
  phe <- subset(phe, select = columns)

################################################################################
#### Format Q and K to proper structure, numeric format
################################################################################

    if (Q_given) {
      Q <- Q[rownames(Q) %in% filtered_accessions,]
      Q <- cbind(Accession = rownames(Q), Q)
      Q <- as.data.frame(Q)
      Q <- Q %>% mutate_at(colnames(Q)[2:length(colnames(Q))], as.numeric)
    }
    if (K_given) {
      K <- K[rownames(K) %in% filtered_accessions,
             colnames(K) %in% filtered_accessions]
      K <- cbind(Accession = rownames(K), K)
      K <- as.data.frame(K)
      K <- K %>% mutate_at(colnames(K)[2:length(colnames(K))], as.numeric)
    }

################################################################################
#### Filter data for the region specified
################################################################################

  if (!is.null(opts$chrom)) {
      map = map[map$V1 == opts$chrom,]
      if (opts$start) {
        map = map[map$V4 >= opts$start,]
      }
      if (opts$end) {
        map = map[map$V4 <= opts$end,]
      }
      markers <- c("IID", names(gt)[(names(gt) %in% map$V2)])
      gt <- subset(gt, select=markers)
    } else {
      if (!is.null(opts$start) | !is.null(opts$end)) {
        stop("Start or end of region specified without specifying chromosome.", call.=FALSE)
      }
    }

################################################################################
#### Tool-specific input validation
################################################################################

    if (tool == "BGLR") {
      if (is.null(opts$chromosome)) {
        stop("Chromosome must be provided.", call.=FALSE)
      }

      if (is.null(opts$start)) {
        stop("Region start must be provided.", call.=FALSE)
      } 

      if (is.null(opts$end)) {
        stop("Region end must be provided.", call.=FALSE)
      }
    }

################################################################################
#### Create outdir and return
################################################################################

    if (!is.null(opts$outdir)) {
      dir.create(file.path(opts$outdir),
                 recursive = T,
                 showWarnings = F)
    }
	  return(list(phe=phe, gt=gt, map=map, K=K, Q=Q))
}