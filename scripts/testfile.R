# 1. Clear environment, garbage collection --------------------------------

rm(list = ls())
gc()


# 2. Load distQTL package, and other packages -----------------------------

library(data.table)
library(distQTL)


# 3. Generate example data table objects, and other inputs ----------------

# Set parameters:
set.seed(31)
nDonors = 100
nSNPs = 50
nGenes = 10
nCells = 1000 * nDonors
p = 8
m = 100
nChromosomes = 1
cisRange = 2e6

# Define cell type groups:
cellTypeGroups = list("B_cells" = c("B"),
                      "Monocytes" = c("Mono"),
                      "CD4T_cells" = c("CD4T"),
                      "CD8T_cells" = c("CD8T"))

# Generate genotype data frame, with donor ID column:
if(TRUE){
  
  probs = rbeta(nGenes, 5, 20)
  genotypeDataTable = data.table::data.table("donorID" = replicate(nDonors, paste(sample(0:9, 2, replace = TRUE),
                                                                                  sample(0:9, 2, replace = TRUE),
                                                                                  sample(0:9, 2, replace = TRUE), sep = "", collapse = "-")),
                                             matrix(rbinom(nDonors * nSNPs, 2, probs), nDonors, nSNPs, byrow = TRUE))
  colnames(genotypeDataTable)[-1] = paste0("rs", replicate(nSNPs, paste(sample(0:9, 7, replace = TRUE), sep = "", collapse = "")))
  
}

# Generate expression data frame, with donor ID column and cell type column:
if(TRUE){
  
  expressionDataTable = data.table::data.table("donorID" = sample(genotypeDataTable$"donorID", nCells, replace = TRUE),
                                               "cellType" = sample(unlist(cellTypeGroups), nCells, replace = TRUE, prob = c(0.2, 0.1, 0.4, 0.3)),
                                               matrix(rpois(nCells * nGenes, 2), nCells, nGenes))
  colnames(expressionDataTable)[3:ncol(expressionDataTable)] = replicate(nGenes, paste(sample(LETTERS, 8, replace = TRUE), sep = "", collapse = ""))
  
}

# Gene info:
if(TRUE){
  
  geneInfo = data.table::data.table("geneID" = colnames(expressionDataTable)[3:ncol(expressionDataTable)],
                                    "chromosome" = sample(1:nChromosomes, nGenes, replace = TRUE),
                                    "start" = round(runif(nGenes, -0.5, 1e6 + 0.5)))
  
}

# SNP info:
if(TRUE){
  
  snpInfo = data.table::data.table("snpID" = colnames(genotypeDataTable)[-1],
                                   "chromosome" = sample(1:nChromosomes, nSNPs, replace = TRUE),
                                   "start" = round(runif(nSNPs, -0.5, 1e6 + 0.5)))
  
}

# Generate covariate data frame, with donor ID column:
if(TRUE){
  
  covariateDataTable = data.table::data.table("donorID" = genotypeDataTable$donor_id,
                                         matrix(rnorm(nDonors * p), nDonors, p))
  
}

# Show all objects:
head(genotypeDataTable[ , 1:5])
head(expressionDataTable[ , 1:5])
head(covariateDataTable)
cellTypeGroups
head(geneInfo)
head(snpInfo)


# 4. Run distQTL ----------------------------------------------------------

t0 = Sys.time()
output = distQTL(genotypeDataTable = genotypeDataTable,
                 expressionDataTable = expressionDataTable,
                 covariateDataTable = covariateDataTable,
                 cellTypeGroups = cellTypeGroups,
                 geneInfo = geneInfo,
                 snpInfo = snpInfo,
                 m = m,
                 cisRange = cisRange)
td = as.numeric(difftime(Sys.time(), t0, units = "secs"))

# Show an example gene from the first cell type group:
output[[1]][2]

# Plot histogram of all p-values:
hist(exp(unlist(output)), breaks = seq(0, 1, 0.005), freq = FALSE,
     main = "distQTL P-value histogram", xlab = "P-value", ylab = "Density",
     border = NA, col = 'gray70')

# Time per model:
sprintf("%.2f milliseconds per model", td / length(unlist(output)) * 1000)


