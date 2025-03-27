#' trimY
#'
#' @param Y 
#' @param lower 
#' @param upper 
#'
#' @returns
#' @export
#'
#' @examples
trimY = function(Y, lower = -Inf, upper = Inf){
  
  wL = which(colSums(Y > lower) == 0)
  wU = which(colSums(Y < upper) == 0)
  w = c(wL, wU)
  Y = if(length(w) == 0) Y else Y[ , -w, drop = FALSE]
  Y

}
