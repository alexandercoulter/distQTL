% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distQTL.R
\name{distQTL}
\alias{distQTL}
\title{distQTL}
\usage{
distQTL(
  genotype = NULL,
  geneExpression = NULL,
  covariates = NULL,
  geneInfo = NULL,
  snpInfo = NULL,
  cellTypeGroups = NULL,
  m = 100,
  cisRange = 1e+05,
  minCells = 10,
  minExpr = 0.01,
  nPermutations = 1
)
}
\arguments{
\item{genotype}{A \code{data.table} object with one row per donor, and \code{1 + nSNPs} columns. One column must be labeled \code{donorID}, which contains unique donor IDs for the data set. The remaining columns must contain donor genotype coding per SNP - \code{0} for homozygous major allele \code{(AA)}, \code{1} for heterozygous genotype \code{(Aa)}, and \code{2} for homozygous minor allele \code{(aa)}. These columns' names should match the values given in the \code{snpID} column in the \code{snpInfo} input object.}

\item{geneExpression}{A \code{data.table} object with one row per cell, and \code{2 + nGenes} columns. One column must be labeled \code{donorID}, which contains the donor ID label per-cell, values matching from the \code{donorID} column. One column must be labeled \code{cellType}, which contains cell type label per-cell. The remaining columns must contain raw gene expression count measurements, e.g. values in \verb{\{0, 1, 2, ...\}}. Sparse vector data types, e.g. from \code{Matrix::sparseVector}, are permitted in these columns. These columns' names should match the values given in the \code{geneID} column in the \code{geneInfo} input object.}

\item{covariates}{A \code{data.table} object with one row per donor, and \code{1 + nCovariates} columns, with precisely as many donors as \code{genotype}; input will be reordered row-wise to match \code{genotype}. One column must be labeled \code{donorID}, which contains unique donor IDs for the data set. The remaining columns should contain desired donor-level covariate vectors, such as demographic information, genotype PCAs, PEER factors, batch information, etc. These columns should contain numeric data only, or data which can be coerced to numeric. Any factor type (or equivalent) data should be expanded into \code{0/1} encoding prior to input. These columns need not have meaningful names.}

\item{geneInfo}{A \code{data.table} object with one row per gene, and \code{3} columns. One column must be labeled \code{geneID}, which contains unique gene identifiers (e.g. gene names, or Ensembl IDs). One column must be labeled \code{chromosome}, which contains unique chromosome identifiers for the genes; the values should match those used in the corresponding \code{chromosome} column of the \code{snpInfo} input object. One column must be labeled \code{start}, which gives the starting base pair locations for the genes in the genome; the reference genome should match the one used in the corresponding \code{start} column in the \code{snpInfo} input object.}

\item{snpInfo}{A \code{data.table} object with one row per SNP, and \code{3} columns. One column must be labeled \code{snpID}, which contains unique SNP identifiers. One column must be labeled \code{chromosome}, which contains unique chromosome identifiers for the SNPs; the values should match those used in the corresponding \code{chromosome} column of the \code{geneInfo} input object. One column must be labeled \code{start}, which gives the starting base pair locations for the SNPs in the genome; the reference genome should match the one used in the corresponding \code{start} column in the \code{geneInfo} input object.}

\item{cellTypeGroups}{An optionally named list object, where each element is a character vector containing cell types desired to be grouped together. These cell types should match from the \code{cellType} column in the \code{geneExpression} input, but need not be exhaustive of those entries. The list entries' names will be used to label cell type group level outputs; we recommend the list entries' names are short and descriptive of the groups they correspond to, e.g. \code{cellTypeGroups = list(B_cells = c("..."), Monocytes = c("..."), ...)}.}

\item{m}{A positive integer (default \code{100}) containing the empirical distribution grid density, i.e. for a given group of cells \code{y} satisfying a donor-cellType-gene combination, the response object is given by \code{quantile(y, seq(from = 0.5 / m, to = 1 - 0.5 / m, length.out = m), type = 1)}.}

\item{cisRange}{A positive integer (default \code{1e5}) specifying the range defining cis-SNPs around each gene's starting location.  cis-SNP range defined by absolute difference in \code{start} values between genes in \code{geneInfo} and SNPs in \code{snpInfo}, conditioned on same \code{chromosome} values.}

\item{minCells}{A positive integer (default \code{10}) specifying the minimum number of cells a donor must have in order to be included in regression.  Evaluated per cell type group.}

\item{minExpr}{A scalar between \eqn{(0, 1]} (default \code{0.01}) giving the minimum required proportion of cells, averaged across donors, that have positive gene expression.  For instance, if \code{minExpr = 0.10}, and there are 2 donors with 7\% and 11\% positive gene expression among their cells, respectively, then the average positive expression across donors is \eqn{(7 + 11)/2 = 9 \leq 10}\%.  This gene would then be considered "low expression", and excluded from distQTL calculations.}

\item{nPermutations}{A non-negative integer (default \code{1}) giving the number of distQTL evaluations per gene-SNP pair, applied under permutation of responses. If nPermutations is positive, \code{distQTL} will also perform p-value correction using these p-values under permutation.}
}
\value{
A nested list structure of distQTL p-values. First list layer splits by cell type group from \code{cellTypeGroups} input; next sub-layer splits by gene from gene columns of \code{geneExpression}. Within each of these sub-layers is a matrix of \eqn{\log_{10}(p)}-values form distQTL. The first column is the raw distQTL \eqn{\log_{10}(p)}-values, individually testing each cis-SNP conditioned on the covariates specified in \code{covariates}. If \code{nPermutations > 0}, the next \code{nPermutations} columns give \eqn{\log_{10}(p)}-values from distQTL applied to data with permuted responses, one column per permutation. Also in that case, the last column gives \eqn{\log_{10}(p)}-values after bulk null-correction.
}
\description{
distQTL
}
\examples{
# See vignette `intro-distQTL`.
}
