#' Wasserstein partial F-test for distribution responses.
#'
#' @param X An $n \times p$ covariate matrix (with no intercept column), whose last column is the one being tested.
#' @param Y An $n \times m$ matrix row-wise containing quantile function response objects.
#' @param lower A numeric scalar (default `-Inf`) setting the lower support box constraint.
#' @param upper A numeric scalar (default `Inf`) setting the upper support box constraint.
#' @param Q0 An optional (i.e. default `NULL`) $n \times m$ matrix containing the fitted quantile functions from the null Fréchet regression model, i.e. the model removing the last column of `X`.
#' @param Qm An optional (i.e. default `NULL`) $n \times m$ matrix containing the marginal Fréchet mean of the response quantile functions; is evaluated as `colMeans(Y)`.
#' @param C_init An optional (i.e. default `NULL`) $n \times (m + 1)$ matrix containing the initial Lagrangian for the QP-problem associated with the Fréchet regression problem. Positive entries correspond to active constraints; zero values correspond to inactive constraints. `NULL` is equivalent to specifying `C_init = matrix(0, n, m + 1)`.
#' @param log.p A boolean value (default `TRUE`) specifying whether the $\log(p)$-value, rather than $p$-value, should be returned from the partial F test.
#'
#' @returns A list with F statistic, p-value ($\log(p)$ if specified in input), and F test degrees of freedom $df_1$ and $df_2$.
#' @export
#'
#' @examples
Wasserstein_F = function(X,
                         Y,
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
  
  if(is.null(Q0)) Q0 = fastfrechet::frechetreg_univar2wass(X = X[ , -p],
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
  
  SYY = Sigma[-p, -p]
  SZY = Sigma[p, -p]
  A = solve(SYY, SZY)
  
  SZ_Y = c(1 - SZY %*% A)
  
  Jt = c(-A, 1)
  C = crossprod(E, c(tcrossprod(Jt, X)^2) * E) / SZ_Y
  
  s2 = mean(C * C)# sum all entries, sum of square eigenvalues
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