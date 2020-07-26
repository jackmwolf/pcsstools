context("calculate_lm is exact")

test_that("calculate_lm betas are exact", {
  ex_data <- cont_data[c("g", "x", "y1")]
  means <- colMeans(ex_data)
  covs <- cov(ex_data)
  n <- nrow(ex_data)

  betas <- calculate_lm(means = means, covs = covs, n = n, add_intercept = TRUE)$beta
  expected_betas <- coef(lm(y1 ~ g + x + 1, data = ex_data))

  expect_equal(betas, expected_betas)
})

test_that("calculate_lm standard errors are exact", {
  ex_data <- cont_data[c("g", "x", "y1")]
  means <- colMeans(ex_data)
  covs <- cov(ex_data)
  n <- nrow(ex_data)

  ses <- calculate_lm(means = means, covs = covs, n = n, add_intercept = TRUE)$sd_beta
  expected_ses <- coef(summary(lm(y1 ~ g + x + 1, data = ex_data)))[, 2]

  expect_equal(ses, expected_ses)
})