#' OneK1K Data Subset
#'
#' These data sets contain a subset of the gene expression and SNP genotype
#' information from the OneK1K PBMC data set.  Included is a random selection
#' of 20 genes and each of their cis-SNPs (i.e. within 200 kpb), along with
#' auxiliary gene and SNP information like labels/names, and chromosome
#' locations.  Each object is a ready-to-use input for the distQTL function.
#'
#' `genotype` and `geneExpression` are data.table objects which contain SNP
#' genotype and library size corrected gene expression, respectively. The first
#' 1/2 columns (again, respectively) contain donor ID; and donor ID and cell
#' type labels.  Other columns have SNP or gene names; see `distQTL`
#' documentation for description of structure, or `intro-distQTL` vignette for
#' header displays of these exact objects.  `covariates` is a data.table object
#' with various demographic and other fixed effects covariates, such as two
#' orthogonal polynomial columns of donor age, binary donor sex labels, and six
#' principal component columns from the full (unprovided) genotype matrix. 
#' `geneInfo` and `snpInfo` are data.table objects which tabulate gene/SNP
#' names/labels, and chromosome locations.  Location entries in these objects
#' are matched and compared to each other for the purpose of identifying 
#' gene-SNP pairs that distQTL model fitting will be performed on; name/label
#' information is matched to the column names of the `geneExpression` and
#' `genotype` objects.  Finally, `cellTypeGroups` is a list object containing
#' the cell type names composing the "B cell" and "Monocyte" groups, which are
#' fit separately by distQTL.
#' 
#' @format Except for `cellTypeGroups`, which is a list, all objects are
#' data.table objects (requiring the `data.table` package).
#' @examples
#' genotype
#' geneExpression
#' covariates
#' geneInfo
#' snpInfo
#' cellTypeGroups
"genotype"
"geneExpression"
"covariates"
"geneInfo"
"snpInfo"
"cellTypeGroups"