#' trimY
#'
#' @param Y A matrix whose rows are monotone non-decreasing and obey `Y >= lower` and `Y <= upper`.
#' @param lower A numeric scalar smaller than `Inf` (default `-Inf`) giving lower distribution support box constraint.  Must be strictly smaller than `upper`.
#' @param upper A numeric scalar larger than `Inf` (default `-Inf`) giving upper distribution support box constraint.  Must be strictly larger than `lower`.
#'
#' @returns A column sub-matrix of `Y`, consisting of those columns not all equal to `lower` and not all equal to `upper`.
#' @export
#'
#' @examples
#' set.seed(31)
#' 
#' # Set number of rows, and quantile function grid density, respectively:
#' n = 5
#' m = 20
#' 
#' # Quantile function grid support:
#' mseq = seq(0.5 / m, 1 - 0.5 / m, length.out = m)
#' 
#' # Generate random quantile functions:
#' Y = t(replicate(n, qpois(mseq, rexp(1, 2))))
#' 
#' # Show Y, which has some columns all zero, the lower bound for Poisson:
#' print(Y)
#' 
#' # Trim the Y matrix:
#' Y_trim = trimY(Y, lower = 0)
#' 
#' # Show trimmed Y, which now has no columns all equal to zero:
#' print(Y_trim)
trimY = function(Y, lower = -Inf, upper = Inf){
  
  wL = which(colSums(Y > lower) == 0)
  wU = which(colSums(Y < upper) == 0)
  w = c(wL, wU)
  Y = if(length(w) == 0) Y else Y[ , -w, drop = FALSE]
  Y

}
