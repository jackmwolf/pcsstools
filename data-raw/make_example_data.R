
invLogit <- function(logit) exp(logit) / (1 + exp(logit))

make_data <- function(response, m, n = 1000, correlated_errors = TRUE) {
  g <- rbinom(n = n, size = 2, prob = 0.3)
  x <- rnorm(n = n)
  betas <- matrix(rnorm(m * 3), nrow = 3)

  fit <- cbind(1, g, x) %*% betas

  sds <- rgamma(n = m, shape = 4, rate = 3)

  if (correlated_errors) {
    require(clusterGeneration)
    Sigma <- clusterGeneration::rcorrmatrix(d = m, alphad = 1)
  } else {
    Sigma <- diag(m)
  }

  # Convert correlation to covariance
  Sigma <- Sigma * (sds %*% t(sds))

  errors <- MASS::mvrnorm(n = n, mu = rep(0, m), Sigma)
  Y <- fit + errors
  if (response == "binary") {
    Y <- matrix(rbinom(n = n * m, size = 1, prob = invLogit(Y)), nrow = n)
  }
  colnames(Y) <- paste0("y", 1:m)

  return(as.data.frame(cbind(g, x, Y)))
}

set.seed(1004)
bin_data  <- make_data(response = "binary", m = 5)
cont_data <- make_data(response = "continuous", m = 5)

save(cont_data, file="data/cont_data.RData")
save(bin_data, file="data/bin_data.RData")


test_data_bin <- make_data(response = "binary", n = 10000, m = 3)
save(test_data_bin, file = "data/test_data_bin.RData")
