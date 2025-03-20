distQTL = function(genotypeMatrix = NULL,
                   expressionMatrix = NULL,
                   covariateMatrix = NULL,
                   cell_cellTypes = NULL,
                   cell_donorIDs = NULL,
                   cellTypeGroups = NULL,
                   m = 100){
  
  # NULL checks:
  if(is.null(genotypeMatrix)) stop("SNP genotype matrix, [n donors] x [n SNPs], must be provided.")
  if(is.null(expressionMatrix)) stop("RNA expression matrix, [n cells] x [n genes], must be provided.")
  if(is.null(covariateMatrix)) stop("Covariate matrix, [n donors] x [n covariates], must be provided.")
  if(is.null(cell_cellTypes)) stop("Vector of cell types per cell, [n cells] long, must be provided.")
  if(is.null(cell_donorIDs)) stop("Vector of donor ids per cell, [n cells] long, must be provided.")
  if(is.null(cellTypeGroups)) print("Cell type labels vector not provided; defaulting to aggregating regardless of cell type.")
  
  mseq = seq(0.5 / m, 1 - 0.5 / m, length.out = m)
  geneLabels = colnames(expressionMatrix)[3:ncol(expressionMatrix)]
  Xcov = as.matrix(covariateMatrix[ , -1])
  
  for(j in 1:length(cellTypeGroups)){
    
    j = 1
    for(i in 1:length(geneLabels)){
      
      i = 1
      Y = expressionMatrix[cellType %in% cellTypeGroups[[j]], as.list(quantile(.SD[[1]], mseq)), by = donorID, .SDcols = geneLabels[i]]
      Y = as.matrix(Y[match(genotypeMatrix$donor_id, donorID), -1])
      
      Q0 = fastfrechet::frechetreg_univar2wass(Xcov, Y, lower = 0, upper = Inf)$Qhat
      Qm = colMeans(Y)
      
      for(k in 1:ncol(genotypeMatrix)){
        
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
          pvalues[[i]][[j]][k] = wass$p_value
          
        }
        
      }
      
    }
    
  }

  return(pvalues)
  
}
