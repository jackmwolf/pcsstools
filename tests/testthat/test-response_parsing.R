context("pcsslm parses the correct response function")

test_that("pcsslm parses linear combinations", {
  ex_data <- pcsstools_example
  
  pcss <- list(
    means = colMeans(ex_data),
    covs = cov(ex_data),
    n = nrow(ex_data)
  )
  
  pcsslm_mod <- pcsslm(y1 + y2 ~ g1 + x1 + x2 + x3, pcss = pcss, phi = c(1, 1))
  model_Z_mod <- model_combo(y1 + y2 ~ g1 + x1 + x2 + x3,
                             phi = c(1,1), n = pcss$n, means = pcss$means,
                             covs = pcss$covs)
  
  # Compare everything but function calls
  expect_equal(pcsslm_mod[-(1)], model_Z_mod[-(1)])
})


test_that("pcsslm parses PCA", {
  ex_data <- pcsstools_example
  
  pcss <- list(
    means = colMeans(ex_data),
    covs = cov(ex_data),
    n = nrow(ex_data)
  )
  
  pcsslm_mod <- pcsslm(y1 + y2 + y3 ~ g1 + x1 + x2 + x3, pcss = pcss, comp = 2)
  model_Z_mod <- model_prcomp(y1 + y2 + y3 ~ g1 + x1 + x2 + x3,
                             comp = 2, n = pcss$n, means = pcss$means,
                             covs = pcss$covs)
  
  # Compare everything but function calls
  expect_equal(pcsslm_mod[-(1)], model_Z_mod[-(1)])
})


test_that("pcsslm parses products", {
  ex_data <- pcsstools_example
  
  pcss <- list(
    means = colMeans(ex_data),
    covs = cov(ex_data),
    n = nrow(ex_data),
    predictors = list(
      g1 = new_predictor_snp(maf = mean(ex_data$g1)/2),
      x1 = new_predictor_normal(mean = mean(ex_data$x1), sd = sd(ex_data$x1))
    )
  )

  
  pcsslm_mod <- pcsslm(y1 * y2 ~ g1 + x1, pcss = pcss)
  model_Z_mod <- model_product(y1 * y2 ~ g1 + x1,
                              n = pcss$n, means = pcss$means,
                              covs = pcss$covs, predictors = pcss$predictors)
  
  # Compare everything but function calls
  expect_equal(pcsslm_mod[-(1)], model_Z_mod[-(1)])
})


test_that("pcsslm parses disjunctions", {
  ex_data <- pcsstools_example
  
  pcss <- list(
    means = colMeans(ex_data),
    covs = cov(ex_data),
    n = nrow(ex_data),
    predictors = list(
      g1 = new_predictor_snp(maf = mean(ex_data$g1)/2),
      x1 = new_predictor_normal(mean = mean(ex_data$x1), sd = sd(ex_data$x1))
    )
  )
  
  
  pcsslm_mod <- pcsslm(y4 | y5 ~ g1 + x1, pcss = pcss)
  model_Z_mod <- model_or(y4 | y5 ~ g1 + x1,
                               n = pcss$n, means = pcss$means,
                               covs = pcss$covs, predictors = pcss$predictors)
  
  # Compare everything but function calls
  expect_equal(pcsslm_mod[-(1)], model_Z_mod[-(1)])
})


test_that("pcsslm parses conjunctions", {
  ex_data <- pcsstools_example
  
  pcss <- list(
    means = colMeans(ex_data),
    covs = cov(ex_data),
    n = nrow(ex_data),
    predictors = list(
      g1 = new_predictor_snp(maf = mean(ex_data$g1)/2),
      x1 = new_predictor_normal(mean = mean(ex_data$x1), sd = sd(ex_data$x1))
    )
  )
  
  
  pcsslm_mod <- pcsslm(y4 & y5 ~ g1 + x1, pcss = pcss)
  model_Z_mod <- model_and(y4 & y5 ~ g1 + x1,
                          n = pcss$n, means = pcss$means,
                          covs = pcss$covs, predictors = pcss$predictors)
  
  # Compare everything but function calls
  expect_equal(pcsslm_mod[-(1)], model_Z_mod[-(1)])
})
