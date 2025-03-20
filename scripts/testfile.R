
# 1. Clear environment, garbage collection --------------------------------

rm(list = ls())
gc()



source("R/distQTL.R")
source("R/Wasserstein_F.R")
Rcpp::sourceCpp("src/scaleX_cpp.cpp")

# Generate example data tables:
set.seed(31)
nDonors = 100
nSNPs = 500
nGenes = 20
nCells = 1000 * nDonors
p = 8
m = 100
nChromosomes = 1
cisRange = 2e6

cellTypeGroups = list("B cells" = c("B"),
                      "Monocytes" = c("Mono"),
                      "CD4T cells" = c("CD4T"),
                      "CD8T cells" = c("CD8T"))

# Generate genotype data frame, with donor ID column:
probs = rbeta(nGenes, 5, 20)
genotypeMatrix = data.table::data.table("donor_id" = replicate(nDonors, paste(sample(0:9, 2, replace = TRUE),
                                                                              sample(0:9, 2, replace = TRUE),
                                                                              sample(0:9, 2, replace = TRUE), sep = "", collapse = "-")),
                                        matrix(rbinom(nDonors * nSNPs, 2, probs), nDonors, nSNPs, byrow = TRUE))
colnames(genotypeMatrix)[-1] = paste0("rs", replicate(nSNPs, paste(sample(0:9, 7, replace = TRUE), sep = "", collapse = "")))

# Generate expression data frame, with donor ID column and cell type column:
expressionMatrix = data.table::data.table("donorID" = sample(genotypeMatrix$"donor_id", nCells, replace = TRUE),
                                          "cellType" = sample(unlist(cellTypeGroups), nCells, replace = TRUE, prob = c(0.2, 0.1, 0.4, 0.3)),
                                          matrix(rpois(nCells * nGenes, 2), nCells, nGenes))
colnames(expressionMatrix)[3:ncol(expressionMatrix)] = replicate(nGenes, paste(sample(LETTERS, 8, replace = TRUE), sep = "", collapse = ""))

# Gene info:
geneInfo = data.table::data.table("geneID" = colnames(expressionMatrix)[3:ncol(expressionMatrix)],
                                  "chromosome" = sample(1:nChromosomes, nGenes, replace = TRUE),
                                  "location" = round(runif(nGenes, -0.5, 1e6 + 0.5)))

# SNP info:
snpInfo = data.table::data.table("snpID" = colnames(genotypeMatrix)[-1],
                                 "chromosome" = sample(1:nChromosomes, nSNPs, replace = TRUE),
                                 "location" = round(runif(nSNPs, -0.5, 1e6 + 0.5)))

# Generate covariate data frame, with donor ID column:
covariateMatrix = data.table::data.table("donorID" = genotypeMatrix$donor_id,
                                         matrix(rnorm(nDonors * p), nDonors, p))

t0 = Sys.time()
distQTL_output = distQTL(genotypeMatrix = genotypeMatrix,
                         expressionMatrix = expressionMatrix,
                         covariateMatrix = covariateMatrix,
                         cellTypeGroups = cellTypeGroups,
                         m = 100,
                         cisRange = cisRange)
difftime(Sys.time(), t0, units = "sec")

hist(exp(unlist(distQTL_output)), breaks = seq(0, 1, 0.001), freq = FALSE)
length(unlist(distQTL_output))
nGenes * length(cellTypeGroups) * nSNPs

distQTL_output$group_1$TYNLKAAG

267 / 40000 * 1.6e6 / 9 / 3600 * 9 * 2