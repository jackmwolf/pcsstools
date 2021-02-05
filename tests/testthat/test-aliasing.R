context("pcsslm aliases coefficients")

test_that("pcsslm aliases coefficients", {
  n <- 100
  y <- rnorm(n)
  x1 <- rnorm(n, 0.5 * y, 3)
  x2 <- 1 + x1
  x3 <- rnorm(n)
  
  dat <- data.frame(y, x1, x2, x3)
  
  pcss <- list(
    means = colMeans(dat),
    covs  = cov(dat),
    n = n
  )
  
  lmpcss <- pcsslm(y ~ x1 + x2 + x3, pcss = pcss)
  lmipd  <- summary(lm(y ~ x1 + x2 + x3, data = dat))
  
  expect_equal(coef(lmpcss), coef(lmipd))
})
