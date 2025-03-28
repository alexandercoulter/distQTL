test_that("Ensure trimY output has no columns equal to lower/upper", {
  set.seed(31)
  
  # Set number of rows, and quantile function grid density, respectively:
  n = 5
  m = 500
  
  # Quantile function grid support:
  mseq = seq(0.5 / m, 1 - 0.5 / m, length.out = m)
  
  # Generate random quantile functions:
  Y = t(replicate(n, qpois(mseq, rexp(1, 1))))
  
  # Constraint Y to between 0 and 3:
  lower = 0
  upper = 3
  Y[Y < lower] = lower
  Y[Y > upper] = upper
  
  # Trim the Y matrix:
  Y_trim = trimY(Y, lower, upper)
  minY = apply(Y_trim, 2, min)
  maxY = apply(Y_trim, 2, max)
  
  # Check that the largest minimum is still less than upper, implying no column
  # has all its values equal to (or larger than) upper:
  expect_lt(max(minY), upper,
            label = "Some column is all 'lower'."
  )
  # Check the smallest maximum is still greater than lower, implying no column
  # has all its values equal to (or smaller than) lower:
  expect_gt(min(maxY), lower,
            label = "Some column is all 'upper'."
  )
  
})
