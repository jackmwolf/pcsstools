
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
bin_data <- make_data(response = "binary", m = 5)
cont_data <- make_data(response = "continuous", m = 5)

save(cont_data, file = "data/cont_data.RData")
save(bin_data, file = "data/bin_data.RData")


test_data_bin <- make_data(response = "binary", n = 10000, m = 3)
save(test_data_bin, file = "data/test_data_bin.RData")






pcsstools_data <- function(n = 1000) {
  Rho_X <- matrix(c(1, 0.5, 0.1, 0.5, 1, 0.25, 0.1,0.25,  1), nrow = 3)
  X <- MASS::mvrnorm(n = n, mu = rep(0, 3), Sigma = Rho_X)
  colnames(X) <- paste0("x", 1:ncol(X))
  
  mafs <- c(0.25, 0.4, 0.15)
  G <- sapply(mafs, function(x) rbinom(n = n, size = 2, prob = x))
  colnames(G) <- paste0("g", 1:ncol(G))
  
  betas <- matrix(c(1, 0.5, 0, 0, 0.5, 0.1,
                    0, -1, -2, 1, 0.5, 0.3,
                    0, 0, 3, 1, 1, 1,
                    0.5, 0.5, -1, 0.1, 0.2, 0.3,
                    1, 0, 1, -0.5, 0.1, -0.5,
                    -0.5, -0.5, -0.2, 0.2, 0.1, 0.3), 
                  nrow = 6, byrow = TRUE)
  
  eta <- cbind(X, G) %*% betas
  Rho_Y <- clusterGeneration::rcorrmatrix(d = 3, alphad = 1)
  sds <- c(rep(0.5, 3))
  Sigma_Y <- Rho_Y * (sds %*% t(sds))
  
  Y_con <- eta[, 1:3] + MASS::mvrnorm(n = n, mu = rep(0, 3), Sigma_Y)
  Y_bin <- matrix(rbinom(n = n * 3, size = 1, prob = invLogit(eta[, 4:6])), nrow = n)
  Y <- cbind(Y_con, Y_bin)
  colnames(Y) <- paste0("y", 1:ncol(Y))
  
  return(as.data.frame(cbind(G, X, Y)))
}

set.seed(1000)
pcsstools_example <- pcsstools_data(1000)
save(pcsstools_example, file = "data/pcsstools_example.RData")
