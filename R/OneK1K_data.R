#' OneK1K genotype data
#' 
#' @docType data
#' @name genotype
#' @usage data(genotype)
#' @format A data.table with 961 rows and 2015 columns:
#' \describe{
#'   \item{donorID}{Vector of anonymized OneK1K donor ID labels.}
#'   \item{other columns}{Donor genotype data (0/1/2) for 2014 different SNPs, with SNP names as column names.}
#' }
#' @author Alexander Coulter \email{coultera@@tamu.edu}
#' @references \url{https://explore.data.humancellatlas.org/projects/f2078d5f-2e7d-4844-8552-f7c41a231e52}
NULL

#' OneK1K gene expression data
#' 
#' @docType data
#' @name geneExpression
#' @usage data(geneExpression)
#' @format A data.table with 144,813 rows and 22 columns:
#' \describe{
#'   \item{donorID}{Vector of anonymized OneK1K donor ID labels, one per cell.}
#'   \item{cellType}{Vector of cell type labels, one per cell.}
#'   \item{other columns}{Donor library size corrected gene expression data for 20 different genes, with gene Ensembl IDs as column names.}
#' }
#' @author Alexander Coulter \email{coultera@@tamu.edu}
#' @references \url{https://explore.data.humancellatlas.org/projects/f2078d5f-2e7d-4844-8552-f7c41a231e52}
NULL

#' OneK1K covariate data
#' 
#' @docType data
#' @name covariates
#' @usage data(covariates)
#' @format A data.table with 961 rows and 10 columns:
#' \describe{
#'   \item{donorID}{Vector of anonymized OneK1K donor ID labels, one per cell.}
#'   \item{age}{Vector of orthogonal polylnomial (first degree) of donor ages.}
#'   \item{age2}{Vector of orthogonal polylnomial (second degree) of donor ages.}
#'   \item{sex}{Vector of binary donor sex class labels (1/2).}
#'   \item{PC_1 through PC_6}{Vectors of first 6 leading genotype principal components.}
#' }
#' @author Alexander Coulter \email{coultera@@tamu.edu}
#' @references \url{https://explore.data.humancellatlas.org/projects/f2078d5f-2e7d-4844-8552-f7c41a231e52}
NULL

#' OneK1K gene information
#' 
#' @docType data
#' @name geneInfo
#' @usage data(geneInfo)
#' @format A data.table with 20 rows and 3 columns:
#' \describe{
#'   \item{geneID}{Vector of gene Ensembl IDs.}
#'   \item{chromosome}{Vector of chromosome IDs (1-22, character).}
#'   \item{start}{Vector of genome sequence start locations (within-chromosome), relative to hg19 reference genome.}
#' }
#' @author Alexander Coulter \email{coultera@@tamu.edu}
#' @references \url{https://explore.data.humancellatlas.org/projects/f2078d5f-2e7d-4844-8552-f7c41a231e52}
NULL

#' OneK1K SNP information
#' 
#' @docType data
#' @name snpInfo
#' @usage data(snpInfo)
#' @format A data.table with 2014 rows and 3 columns:
#' \describe{
#'   \item{snpID}{Vector of SNP names as provided from OneK1K data pull.}
#'   \item{chromosome}{Vector of chromosome IDs (1-22, character).}
#'   \item{start}{Vector of genome sequence start locations (within-chromosome), relative to hg19 reference genome.}
#' }
#' @author Alexander Coulter \email{coultera@@tamu.edu}
#' @references \url{https://explore.data.humancellatlas.org/projects/f2078d5f-2e7d-4844-8552-f7c41a231e52}
NULL

#' OneK1K cell type groups
#' 
#' @docType data
#' @name cellTypeGroups
#' @usage data(cellTypeGroups)
#' @format A named list with two entries:
#' \describe{
#'   \item{B}{Vector of B cell type names ("naive B cell", "memory B cell", and "transitional stage B cell").}
#'   \item{Mono}{Vector of monocyte type names ("CD14-positive monocyte", and "CD14-low, CD16-positive monocyte").}
#' }
#' @author Alexander Coulter \email{coultera@@tamu.edu}
#' @references \url{https://explore.data.humancellatlas.org/projects/f2078d5f-2e7d-4844-8552-f7c41a231e52}
NULL
