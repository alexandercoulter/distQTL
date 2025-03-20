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