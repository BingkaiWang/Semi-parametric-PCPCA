# simulation 1: testing the validity of the error measure under different asymptotics
# environment setup ------------
set.seed(123)
library(multigroup) # running Flury's algorithm
source("simulations/pcpc_test.R")  # loading required functions and packages 
p <-20
k <- 10
lambda <- exp(seq(0, 10, by = 0.5))[1:20]
n_sim <- 1000
print("simulation 1: testing the validity of the error measure under different asymptotics")
# parallelism
library(foreach)
library(doParallel)
cl <- makeCluster(32)
registerDoParallel(cl)

# asymptotics of n ----------
ns <- 50
n_pool <- c(50, 100, 500, 1000)
lambda_common <- sample(lambda, k, replace = F)
lambda_uncommon <- setdiff(lambda, lambda_common)
gamma_true <- svd(matrix(rnorm(p^2), p, p))$u[,1:k]
sim1_n <- foreach(i=1:n_sim, .combine = cbind, .packages=c('tidyverse','pracma', 'multigroup', 'MASS')) %dopar% {
  err_k <- rep(NA, 4)
  err_kplus1 <- rep(NA, 4)
  for(iter in 1:4){
    n <- n_pool[iter]
    sim_cov <- map(1:n, function(j){
      temp <- gramSchmidt(cbind(gamma_true, matrix(rnorm(p*(p-k)), p, p-k)))$Q[,(k+1):p]
      gamma_true %*% diag(map_dbl(lambda_common, ~(.)), nrow = k) %*% t(gamma_true) +
        temp %*% diag(map_dbl(lambda_uncommon, ~(.)), nrow = p-k) %*% t(temp)
    })
    sim_data <- map(sim_cov, ~cov(MASS::mvrnorm(ns, rep(0,p), .)))
    avg_sim_data <- avg_mat(sim_data)
    pc_semi <- svd(avg_sim_data)$u
    recon_sim_data_semi <- map(sim_data, ~t(pc_semi) %*% . %*% pc_semi)
    eig_avg_sim_data <- eig(avg_sim_data)
    err_semi <- map_dbl(1:p, function(i){
      mean(map_dbl(recon_sim_data_semi, ~sum(.[i,-i]^2/(eig_avg_sim_data[i]*eig_avg_sim_data[-i])))) / (p-1)
    })
    err_k[iter] <- sort(err_semi)[k]
    err_kplus1[iter] <- sort(err_semi)[k+1]
  }
  c(err_k, err_kplus1)
}

n <- 10000
sim_cov <- map(1:n, function(j){
  temp <- gramSchmidt(cbind(gamma_true, matrix(rnorm(p*(p-k)), p, p-k)))$Q[,(k+1):p]
  gamma_true %*% diag(map_dbl(lambda_common, ~(.)), nrow = k) %*% t(gamma_true) +
    temp %*% diag(map_dbl(lambda_uncommon, ~(.)), nrow = p-k) %*% t(temp)
})
sim_data <- map(sim_cov, ~cov(MASS::mvrnorm(ns, rep(0,p), .)))
avg_sim_data <- avg_mat(sim_data)
pc_semi <- svd(avg_sim_data)$u
recon_sim_data_semi <- map(sim_data, ~t(pc_semi) %*% . %*% pc_semi)
eig_avg_sim_data <- eig(avg_sim_data)
err_semi <- map_dbl(1:p, function(i){
  mean(map_dbl(recon_sim_data_semi, ~sum(.[i,-i]^2/(eig_avg_sim_data[i]*eig_avg_sim_data[-i])))) / (p-1)
})
sort(err_semi)[k+1]

# asymptotics of ns ----------
n <- 50
ns_pool <- c(50, 100, 500, 1000)
lambda_common <- sample(lambda, k, replace = F)
lambda_uncommon <- setdiff(lambda, lambda_common)
gamma_true <- svd(matrix(rnorm(p^2), p, p))$u[,1:k]
sim_cov <- map(1:n, function(j){
  temp <- gramSchmidt(cbind(gamma_true, matrix(rnorm(p*(p-k)), p, p-k)))$Q[,(k+1):p]
  gamma_true %*% diag(map_dbl(lambda_common, ~(.)), nrow = k) %*% t(gamma_true) +
    temp %*% diag(map_dbl(lambda_uncommon, ~(.)), nrow = p-k) %*% t(temp)
})
sim1_ns <- foreach(i=1:n_sim, .combine = cbind, .packages=c('tidyverse','pracma', 'multigroup', 'MASS')) %dopar% {
  err_k <- rep(NA, 4)
  err_kplus1 <- rep(NA, 4)
  for(iter in 1:4){
    ns <- ns_pool[iter]
    sim_data <- map(sim_cov, ~cov(MASS::mvrnorm(ns, rep(0,p), .)))
    avg_sim_data <- avg_mat(sim_data)
    pc_semi <- svd(avg_sim_data)$u
    recon_sim_data_semi <- map(sim_data, ~t(pc_semi) %*% . %*% pc_semi)
    eig_avg_sim_data <- eig(avg_sim_data)
    err_semi <- map_dbl(1:p, function(i){
      mean(map_dbl(recon_sim_data_semi, ~sum(.[i,-i]^2/(eig_avg_sim_data[i]*eig_avg_sim_data[-i])))) / (p-1)
    })
    err_k[iter] <- sort(err_semi)[k]
    err_kplus1[iter] <- sort(err_semi)[k+1]
  }
  c(err_k, err_kplus1)
}

saveRDS(list(sim1_n = sim1_n, sim1_ns = sim1_ns), "simulation1.rds")

# fiding the theoretical limit of non-CPC estimate
avg_sim_data <- avg_mat(sim_cov)
pc_semi <- svd(avg_sim_data)$u
recon_sim_data_semi <- map(sim_cov, ~t(pc_semi) %*% . %*% pc_semi)
eig_avg_sim_data <- eig(avg_sim_data)
err_semi <- map_dbl(1:p, function(i){
  mean(map_dbl(recon_sim_data_semi, ~sum(.[i,-i]^2/(eig_avg_sim_data[i]*eig_avg_sim_data[-i])))) / (p-1)
})
sort(err_semi)[k+1]
