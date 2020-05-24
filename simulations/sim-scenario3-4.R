# simulation 3-4: method comparison given k under p = 100, non-Gaussian distribution and
# assuming CPC are associated with largest eigenvalues.
# environment setup ------------
set.seed(123)
library(multigroup) # running Flury's algorithm
source("R/pcpc_test.R")  # loading required functions and packages 
p <-100
k <- 20
lambda <- exp(seq(0, 10, by = 0.1))[1:100]
lambda_uncommon <- lambda[1:80]
lambda_common <- lambda[81:100]
n <- 1000
ns <- 1000
n_sim <- 1000
print("simulation 3-4: method comparison given k under p = 100, Gaussian distribution and fixed CPC rankings")
# parallelism
library(foreach)
library(doParallel)
cl <- makeCluster(32)
registerDoParallel(cl)

sim3.4 <- foreach(i=1:n_sim, .combine = cbind, .packages=c('tidyverse','pracma', 'multigroup', 'MASS')) %dopar% {
  # generate data
  gamma_true <- svd(matrix(rnorm(p^2), p, p))$u[,1:k]
  sim_cov <- map(1:n, function(j){
    temp <- gramSchmidt(cbind(gamma_true, matrix(rnorm(p*(p-k)), p, p-k)))$Q[,(k+1):p]
    gamma_true %*% diag(map_dbl(lambda_common, ~(.)), nrow = k) %*% t(gamma_true) +
      temp %*% diag(map_dbl(lambda_uncommon, ~(.)), nrow = p-k) %*% t(temp)
  })
  sim_data <- map(sim_cov, function(j){
    temp <- matrix(rgamma(ns * p, shape = 0.04, rate = 0.2), nrow = ns) - 0.2
    chol_up <- chol(j)
    temp %*% chol_up
  })
  
  # semi-parametric
  semi_data <- map(sim_data, ~cov(.))
  avg_semi_data <- avg_mat(semi_data)
  pc_semi <- svd(avg_semi_data)$u
  recon_semi_data <- map(semi_data, ~t(pc_semi) %*% . %*% pc_semi)
  eig_avg_semi_data <- eig(avg_semi_data)
  err_semi <- map_dbl(1:p, function(i){
    mean(map_dbl(recon_semi_data, ~sum(.[i,-i]^2/(eig_avg_semi_data[i]*eig_avg_semi_data[-i])))) / (p-1)
  })
  pc_hat_semi <- pc_semi[,order(err_semi)[1:k]]

  # semi-parametric assuming k is unknown
  k_hat <- length(pcpc_test(pc_semi, semi_data, ns)$cpc_index)
  pc_hat_semi_k_unknown <- pc_semi[,order(err_semi)[1:k_hat]]
  pc_hat_semi_k_unknown <- rep(NA, p)
  
  # Flury's
  sim_data_stack <- NULL
  sim_data_group <- NULL
  for(j in 1:n){
    sim_data_stack <- rbind(sim_data_stack, sim_data[[j]])
    sim_data_group <- c(sim_data_group, rep(j, ns))
  }
  pc_flury <- FCPCA(Data = sim_data_stack, Group = sim_data_group, graph = F)$loadings.common
  recon_flury_data <- map(sim_data, ~ t(pc_flury) %*% cov(.) %*% pc_flury)
  var_explained_flury <- map_dbl(1:p, function(i){
    mean(map_dbl(recon_flury_data, ~.[i,i]))
  })
  pc_hat_flury <- pc_flury[,order(var_explained_flury, decreasing = T)[1:k]]
  pc_hat_flury <- rep(NA, p)
  
  # PVD
  subject_pc <- map(sim_data, ~(svd(.)$v)[,1:k])
  subject_pc_stack <- NULL
  for(j in 1:n){
    subject_pc_stack <- cbind(subject_pc_stack, subject_pc[[j]])
  }
  pc_hat_pvd <- svd(subject_pc_stack)$u[,1:k]

  # result
  c(mean(apply(abs(t(gamma_true) %*% pc_hat_semi),1,max)),
    mean(apply(abs(t(gamma_true) %*% pc_hat_semi_k_unknown),1,max)),
    mean(apply(abs(t(gamma_true) %*% pc_hat_flury),1,max)),
    mean(apply(abs(t(gamma_true) %*% pc_hat_pvd),1,max)))
}

saveRDS(sim3.4, "simulation3-4.rds")