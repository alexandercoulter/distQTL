#' distQTL
#'
#' @param genotype A `data.table` object with one row per donor, and `1 + nSNPs` columns. One column must be labeled `donorID`, which contains unique donor IDs for the data set. The remaining columns must contain donor genotype coding per SNP - `0` for homozygous major allele `(AA)`, `1` for heterozygous genotype `(Aa)`, and `2` for homozygous minor allele `(aa)`. These columns' names should match the values given in the `snpID` column in the `snpInfo` input object.
#' @param geneExpression A `data.table` object with one row per cell, and `2 + nGenes` columns. One column must be labeled `donorID`, which contains the donor ID label per-cell, values matching from the `donorID` column. One column must be labeled `cellType`, which contains cell type label per-cell. The remaining columns must contain raw gene expression count measurements, e.g. values in `{0, 1, 2, ...}`. Sparse vector data types, e.g. from `Matrix::sparseVector`, are permitted in these columns. These columns' names should match the values given in the `geneID` column in the `geneInfo` input object.
#' @param covariates A `data.table` object with one row per donor, and `1 + nCovariates` columns, with precisely as many donors as `genotype`; input will be reordered row-wise to match `genotype`. One column must be labeled `donorID`, which contains unique donor IDs for the data set. The remaining columns should contain desired donor-level covariate vectors, such as demographic information, genotype PCAs, PEER factors, batch information, etc. These columns should contain numeric data only, or data which can be coerced to numeric. Any factor type (or equivalent) data should be expanded into `0/1` encoding prior to input. These columns need not have meaningful names.
#' @param geneInfo A `data.table` object with one row per gene, and `3` columns. One column must be labeled `geneID`, which contains unique gene identifiers (e.g. gene names, or Ensembl IDs). One column must be labeled `chromosome`, which contains unique chromosome identifiers for the genes; the values should match those used in the corresponding `chromosome` column of the `snpInfo` input object. One column must be labeled `start`, which gives the starting base pair locations for the genes in the genome; the reference genome should match the one used in the corresponding `start` column in the `snpInfo` input object.
#' @param snpInfo A `data.table` object with one row per SNP, and `3` columns. One column must be labeled `snpID`, which contains unique SNP identifiers. One column must be labeled `chromosome`, which contains unique chromosome identifiers for the SNPs; the values should match those used in the corresponding `chromosome` column of the `geneInfo` input object. One column must be labeled `start`, which gives the starting base pair locations for the SNPs in the genome; the reference genome should match the one used in the corresponding `start` column in the `geneInfo` input object.
#' @param cellTypeGroups An optionally named list object, where each element is a character vector containing cell types desired to be grouped together. These cell types should match from the `cellType` column in the `geneExpression` input, but need not be exhaustive of those entries. The list entries' names will be used to label cell type group level outputs; we recommend the list entries' names are short and descriptive of the groups they correspond to, e.g. `cellTypeGroups = list(B_cells = c("..."), Monocytes = c("..."), ...)`.
#' @param m A positive integer (default `100`) containing the empirical distribution grid density, i.e. for a given group of cells `y` satisfying a donor-cellType-gene combination, the response object is given by `quantile(y, seq(from = 0.5 / m, to = 1 - 0.5 / m, length.out = m), type = 1)`.
#' @param cisRange A positive integer (default `1e5`) specifying the range defining cis-SNPs around each gene's starting location.  cis-SNP range defined by absolute difference in `start` values between genes in `geneInfo` and SNPs in `snpInfo`, conditioned on same `chromosome` values.
#' @param minCells A positive integer (default `10`) specifying the minimum number of cells a donor must have in order to be included in regression.  Evaluated per cell type group.
#' @param minExpr A scalar between \eqn{(0, 1]} (default `0.01`) giving the minimum required proportion of cells, averaged across donors, that have positive gene expression.  For instance, if `minExpr = 0.10`, and there are 2 donors with 7% and 11% positive gene expression among their cells, respectively, then the average positive expression across donors is \eqn{(7 + 11)/2 = 9 \leq 10}%.  This gene would then be considered "low expression", and excluded from distQTL calculations.
#'
#' @returns A nested list structure of distQTL p-values. First list layer splits by cell type group from `cellTypeGroups` input; next sub-layer splits by gene from gene columns of `geneExpression`; within each of these sub-layers is a vector of \eqn{\log_{10}(p)} values form distQTL, individually testing each cis-SNP conditioned on the covariates specified in `covariates`.
#' @export
#' 
#' @importFrom stats quantile
#' @import data.table
#' @import utils
#'
#' @examples
#' # See vignette `intro-distQTL`.
#' @import data.table
#' @import fastfrechet
distQTL = function(genotype = NULL,
                   geneExpression = NULL,
                   covariates = NULL,
                   geneInfo = NULL,
                   snpInfo = NULL,
                   cellTypeGroups = NULL,
                   m = 100,
                   cisRange = 1e5,
                   minCells = 10,
                   minExpr = 0.01){
  
  # Create "fake" global variables to trick parser that cannot properly
  # interpret data.table syntax:
  cellType <- NULL
  donorID <- NULL
  
  # NULL checks:
  if(is.null(genotype)) stop("SNP genotype data.table must be provided.")
  if(is.null(geneExpression)) stop("RNA expression data.table, [n cells] x [n genes], must be provided.")
  if(is.null(covariates)) stop("Covariate matrix, [n donors] x [n covariates], must be provided.")
  if(is.null(cellTypeGroups)) print("Cell type labels vector not provided; defaulting to aggregating regardless of cell type.")
  
  # Extract dimensions:
  nGenes = ncol(geneExpression) - 2
  nSNPs = ncol(genotype) - 1
  nGroups = length(cellTypeGroups)
  
  # Set quantile function support grid in (0, 1):
  mseq = seq(0.5 / m, 1 - 0.5 / m, length.out = m)
  
  # Extract gene and SNP names from geneExpression and genotype, respectively:
  geneNames = colnames(geneExpression)[3:ncol(geneExpression)]
  snpNames = colnames(genotype)[-1]
  
  # Extract covariate information:
  Xcov = as.matrix(covariates[ , -1])
  XdonorID = covariates$donorID
  
  # Match genotype row-wise to covariates:
  genotype = genotype[match(XdonorID, donorID), ]
  
  # Initialize p-value list object, and name the elements:
  pvalues = vector(mode = "list", length = length(cellTypeGroups))
  names(pvalues) = if(!is.null(names(cellTypeGroups))) names(cellTypeGroups) else paste("group", seq_len(length(cellTypeGroups)), sep = "_")
  
  for(j in seq_len(nGroups)){
    
    # j = 2
    
    # Filter out low-expression genes:
    numerator = geneExpression[cellType %in% cellTypeGroups[[j]], lapply(.SD, function(x) if(length(x) == 0) 0 else mean(x > 0)), by = donorID, .SDcols = geneNames]
    denominator = geneExpression[cellType %in% cellTypeGroups[[j]], lapply(.SD, function(x) length(x) >= minCells), by = donorID, .SDcols = geneNames]
    
    # Extract only those donors that have cells in this cell type group:
    Xcov_j = Xcov[XdonorID %in% numerator$donorID, ]
    XdonorID_j = XdonorID[XdonorID %in% numerator$donorID]
    
    # Order with respect to XdonorID_j:
    numerator = numerator[match(XdonorID_j, donorID), ]
    denominator = denominator[match(XdonorID_j, donorID), ]
    
    lowExpr = (colSums(numerator[ , -1] * denominator[ , -1]) / colSums(denominator[ , -1])) < minExpr
    keepGenes = geneNames[which(!lowExpr)]
    
    # Initialize the j^th cell type group p-value list, and name the elements
    # with gene names:
    pvalues[[j]] = vector(mode = "list", length = length(keepGenes))
    names(pvalues[[j]]) = keepGenes
    
    for(i in seq_len(length(keepGenes))){
      
      # i = 1
      # Find cis-SNPs:
      wGI = which(geneInfo$geneID == keepGenes[i])
      cisSNP = (snpInfo$chromosome == geneInfo$chromosome[wGI]) & (abs(geneInfo$start[wGI] - snpInfo$start) <= cisRange)
      keepSNPs = snpNames[which(cisSNP)]
      
      # Initialize p-value storage:
      pvalues[[j]][[i]] = rep(NA, length(keepSNPs))
      names(pvalues[[j]][[i]]) = keepSNPs
      
      if(length(keepSNPs) > 0){
        
        # Extract gene expression information:
        Y = geneExpression[cellType %in% cellTypeGroups[[j]], as.list(quantile(.SD[[1]], mseq)), by = donorID, .SDcols = keepGenes[i]]
        Y = as.matrix(Y[match(XdonorID_j, donorID), -1])
        
        # Remove rows of Xcov and Y corresponding to donors with low cell count:
        wEnoughCells = which(denominator[[which(colnames(denominator) == keepGenes[i])]])
        Xcov_i = Xcov_j[wEnoughCells, , drop = FALSE]
        Y = Y[wEnoughCells, , drop = FALSE]
        
        # Trim columns from Y which are all equal to lower or upper:
        Y = trimY(Y, lower = 0, upper = Inf)
        
        # Fit null model quantile functions and central quantile function:
        Q0 = fastfrechet::frechetreg_univar2wass(Xcov_i, Y, lower = 0, upper = Inf)$Qhat
        Qm = colMeans(Y)
        
        # Loop over cis-SNPs:
        for(k in seq_len(length(keepSNPs))){
          
          # Covariate matrix:
          X = cbind(Xcov_i, genotype[[which(colnames(genotype) == keepSNPs[k])]][wEnoughCells])
          
          # Run Wasserstein F test:
          wass = Wasserstein_F(X = X,
                               Y = Y,
                               lower = 0,
                               upper = Inf,
                               Q0 = Q0,
                               Qm = Qm,
                               C_init = NULL,
                               log.p = TRUE)
          pvalues[[j]][[i]][k] = wass$p_value / log(10)
          
        }
        
      }
      
    }
    
  }
  
  return(pvalues)
  
}
