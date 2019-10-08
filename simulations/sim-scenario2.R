# simulation 2: testing the finding-k algorithm
# environment setup ------------
set.seed(123)
library(multigroup) # running Flury's algorithm
source("simulations/pcpc_test.R")  # loading required functions and packages 
p <-20
lambda <- exp(seq(0, 10, by = 0.5))[1:20]
n <- 100
ns <- 100
n_sim <- 1000
print("simulation 2: testing the finding-k algorithm")
# parallelism
library(foreach)
library(doParallel)
cl <- makeCluster(32)
registerDoParallel(cl)

# Gaussian distirbution ----------
k_pool <- c(1, 10, 19)
sim2.1 <- foreach(i=1:n_sim, .combine = cbind, .packages=c('tidyverse','pracma', 'multigroup', 'MASS')) %dopar% {
  typeIerror <- rep(NA, 3)
  power <- rep(NA, 3)
  for(iter in 1:3){
    k <- k_pool[iter]
    lambda_common <- sample(lambda, k, replace = F)
    lambda_uncommon <- setdiff(lambda, lambda_common)
    gamma_true <- svd(matrix(rnorm(p^2), p, p))$u[,1:k]
    sim_cov <- map(1:n, function(j){
      temp <- gramSchmidt(cbind(gamma_true, matrix(rnorm(p*(p-k)), p, p-k)))$Q[,(k+1):p]
      gamma_true %*% diag(map_dbl(lambda_common, ~(.)), nrow = k) %*% t(gamma_true) +
        temp %*% diag(map_dbl(lambda_uncommon, ~(.)), nrow = p-k) %*% t(temp)
    })
    sim_data <- map(sim_cov, ~cov(MASS::mvrnorm(ns, rep(0,p), .)))
    avg_sim_data <- avg_mat(sim_data)
    pc_semi <- svd(avg_sim_data)$u
    k_hat <- length(pcpc_test(pc_semi, sim_data, ns)$cpc_index)
    typeIerror[iter] <- (k_hat > k)
    power[iter] <- (k_hat == k)
  }
  c(typeIerror, power)
}
saveRDS(sim2.1, "simulation2-1.rds")


# non-Gaussian distribution ------
k_pool <- c(1, 10, 19)
sim2.2 <- foreach(i=1:n_sim, .combine = cbind, .packages=c('tidyverse','pracma', 'multigroup', 'MASS')) %dopar% {
  typeIerror <- rep(NA, 3)
  power <- rep(NA, 3)
  for(iter in 1:3){
    k <- k_pool[iter]
    lambda_common <- sample(lambda, k, replace = F)
    lambda_uncommon <- setdiff(lambda, lambda_common)
    gamma_true <- svd(matrix(rnorm(p^2), p, p))$u[,1:k]
    sim_cov <- map(1:n, function(j){
      temp <- gramSchmidt(cbind(gamma_true, matrix(rnorm(p*(p-k)), p, p-k)))$Q[,(k+1):p]
      gamma_true %*% diag(map_dbl(lambda_common, ~(.)), nrow = k) %*% t(gamma_true) +
        temp %*% diag(map_dbl(lambda_uncommon, ~(.)), nrow = p-k) %*% t(temp)
    })
    sim_data <- map(sim_cov, function(j){
      temp <- matrix(rgamma(ns * p, shape = 0.04, rate = 0.2), nrow = ns) - 0.2
      chol_up <- chol(j)
      cov(temp %*% chol_up)
    })
    avg_sim_data <- avg_mat(sim_data)
    pc_semi <- svd(avg_sim_data)$u
    k_hat <- length(pcpc_test(pc_semi, sim_data, ns)$cpc_index)
    typeIerror[iter] <- (k_hat > k)
    power[iter] <- (k_hat == k)
  }
  c(typeIerror, power)
}
saveRDS(sim2.2, "simulation2-2.rds")
