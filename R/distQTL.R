distQTL = function(genotypeMatrix = NULL,
                   expressionMatrix = NULL,
                   covariateMatrix = NULL,
                   cellTypeGroups = NULL,
                   m = 100,
                   cisRange = 1e5){
  
  # NULL checks:
  if(is.null(genotypeMatrix)) stop("SNP genotype matrix, [n donors] x [n SNPs], must be provided.")
  if(is.null(expressionMatrix)) stop("RNA expression matrix, [n cells] x [n genes], must be provided.")
  if(is.null(covariateMatrix)) stop("Covariate matrix, [n donors] x [n covariates], must be provided.")
  if(is.null(cellTypeGroups)) print("Cell type labels vector not provided; defaulting to aggregating regardless of cell type.")
  
  # Extract dimensions:
  nGenes = ncol(expressionMatrix) - 2
  nSNPs = ncol(genotypeMatrix) - 1
  nGroups = length(cellTypeGroups)
  
  # Set quantile function support grid in (0, 1):
  mseq = seq(0.5 / m, 1 - 0.5 / m, length.out = m)
  
  # Extract gene and SNP names from expression and genotype matrix, respectively:
  geneNames = colnames(expressionMatrix)[3:ncol(expressionMatrix)]
  snpNames = colnames(genotypeMatrix)[-1]
  
  # Extract covariate information:
  Xcov = as.matrix(covariateMatrix[ , -1])
  
  # Initialize p-value list object, and name the elements:
  pvalues = vector(mode = "list", length = length(cellTypeGroups))
  names(pvalues) = paste("group", seq_len(length(cellTypeGroups)), sep = "_")
  
  for(j in seq_len(nGroups)){
    
    # j = 1
    
    # Initialize the j^th cell type group p-value list, and name the elements
    # with gene names:
    pvalues[[j]] = vector(mode = "list", length = length(geneNames))
    names(pvalues[[j]]) = geneNames
    
    for(i in seq_len(nGenes)){
      
      # i = 1
      # Extract gene expression information:
      Y = expressionMatrix[cellType %in% cellTypeGroups[[j]], as.list(quantile(.SD[[1]], mseq)), by = donorID, .SDcols = geneNames[i]]
      Y = as.matrix(Y[match(genotypeMatrix$donor_id, donorID), -1])
      
      # Fit null model quantile functions and central quantile function:
      Q0 = fastfrechet::frechetreg_univar2wass(Xcov, Y, lower = 0, upper = Inf)$Qhat
      Qm = colMeans(Y)
      
      pvalues[[j]][[i]] = rep(NA, ncol(genotypeMatrix) - 1)
      names(pvalues[[j]][[i]]) = snpNames
      
      for(k in seq_len(nSNPs)){
        
        # k = 2
        # Check the SNP is a cis-SNP for the given gene:
        if(geneInfo$chromosome[i] == snpInfo$chromosome[k] && (abs(geneInfo$location[i] - snpInfo$location[k]) <= cisRange)){
          
          # Covariate matrix:
          X = cbind(Xcov, genotypeMatrix[[k + 1]])
          
          # Run Wasserstein F test:
          wass = Wasserstein_F(X = X,
                               Y = Y,
                               lower = 0,
                               upper = Inf,
                               Q0 = Q0,
                               Qm = Qm,
                               C_init = NULL,
                               log.p = TRUE)
          pvalues[[j]][[i]][k] = wass$p_value
          
        }
        
      }
      
    }
    
  }

  return(pvalues)
  
}
