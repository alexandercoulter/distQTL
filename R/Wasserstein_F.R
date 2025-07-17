#' Wasserstein partial F-test for distribution responses.
#'
#' @param X An \eqn{n \times p} covariate matrix (with no intercept column), whose last column is the one being tested.
#' @param Y An \eqn{n \times m} matrix row-wise containing quantile function response objects.
#' @param test A numeric integer (default `NULL`) between `1` and `ncol(X)` indicating which covariate to test; `NULL` means the last covariate of `X` is tested.
#' @param lower A numeric scalar (default `-Inf`) setting the lower support box constraint.
#' @param upper A numeric scalar (default `Inf`) setting the upper support box constraint.
#' @param Q0 An optional (i.e. default `NULL`) \eqn{n \times m} matrix containing the fitted quantile functions from the null Fréchet regression model, i.e. the model removing the last column of `X`.
#' @param Qm An optional (i.e. default `NULL`) \eqn{n \times m} matrix containing the marginal Fréchet mean of the response quantile functions; is evaluated as `colMeans(Y)`.
#' @param C_init An optional (i.e. default `NULL`) \eqn{n \times (m + 1)} matrix containing the initial Lagrangian for the QP-problem associated with the Fréchet regression problem. Positive entries correspond to active constraints; zero values correspond to inactive constraints. `NULL` is equivalent to specifying `C_init = matrix(0, n, m + 1)`.
#' @param log.p A boolean value (default `TRUE`) specifying whether the \eqn{\log(p)}-value, rather than \eqn{p}-value, should be returned from the partial F test.
#'
#' @returns A list with F statistic, p-value (\eqn{\log(p)} if specified in input), and F test degrees of freedom \eqn{df_1} and \eqn{df_2}.
#' @export
#' 
#' @importFrom stats pf
#'
#' @examples
#' # Generate `zinbinom` distributions:
#' n <- 100 # number of samples - nrow(X) and nrow(Y).
#' p <- 4   # number of covariates - ncol(X).
#' m <- 100 # EQF grid density - ncol(Y).
#' mseq <- seq(1 / (2 * m), 1 - 1 / (2 * m), length.out = m)
#'
#' set.seed(31)
#' mydata <- fastfrechet::generate_zinbinom_qf(n = n, p = p, m = m)
#'
#' X <- mydata$X # (n x p) matrix of covariates
#' Y <- mydata$Y # (n x m) matrix of EQFs, stored row-wise
#' 
#' # Run Wasserstein F test on last covariate:
#' output <- Wasserstein_F(X, Y, lower = 0, log.p = FALSE)
#' 
#' # p-value is low:
#' output$p_value
Wasserstein_F = function(X,
                         Y,
                         test = NULL,
                         lower = -Inf,
                         upper = Inf,
                         Q0 = NULL,
                         Qm = NULL,
                         C_init = NULL,
                         log.p = TRUE){
  
  # Get dimensions:
  n = nrow(X)
  m = ncol(Y)
  p = ncol(X)
  
  if(is.null(test)) test = p
  if((test %% 1) != 0 || test < 1 || test > p) stop("test must be NULL or an integer between 1 and p, inclusive.")
  
  if(is.null(Q0)) Q0 = fastfrechet::frechetreg_univar2wass(X = X[ , -test],
                                                           Y = Y,
                                                           C_init = C_init,
                                                           lower = lower,
                                                           upper = upper)$Qhat
  if(is.null(Qm)) Qm = colMeans(Y)
  Qa = fastfrechet::frechetreg_univar2wass(X = X,
                                           Y = Y,
                                           C_init = C_init,
                                           lower = lower,
                                           upper = upper)$Qhat
  
  # Get test statistic:
  RR0 = sum(rowMeans((Q0 - Y)^2))
  RR1 = sum(rowMeans((Qa - Y)^2))
  Fstat = (RR0 - RR1) / (RR1 / (n - p))
  
  # Get E0:
  E = Qa - Y
  
  # Center and scale inputs:
  X = scaleX_cpp(X)
  Sigma = crossprod(X)
  
  SYY = Sigma[-test, -test]
  SZY = Sigma[test, -test]
  A = solve(SYY, SZY)
  
  SZ_Y = c(n - SZY %*% A)
  
  Jt = c(-A, 1)
  C = crossprod(E, c(tcrossprod(Jt, X)^2) * E) / SZ_Y
  
  s2 = mean(C * C) # sum all entries, sum of square eigenvalues
  s1 = mean(diag(C)) # trace, sum of eigenvalues
  
  a = s2 / s1
  b = s1 / a
  
  f1 = b
  f2 = b * (n - p)
  
  pval = if(log.p) pf(Fstat, f1, f2, lower.tail = F, log.p = TRUE) else pf(Fstat, f1, f2, lower.tail = F, log.p = FALSE)
  
  return(list('Fstat' = Fstat,
              'p_value' = pval,
              'df1' = f1,
              'df2' = f2))
  
}