library(tidyverse)
library(pracma) # eigenvalue of a matrix; gramschimdt
library(MASS)
library(Matrix)

pcpc_test <- function(gamma, cov_data, ns, m = 100){
  reconstruct_cov_data <- map(cov_data, ~t(gamma) %*% . %*% gamma)
  p <- ncol(gamma)
  n <- length(cov_data)
  avg_cov_data <- avg_mat(cov_data)
  eig_avg_cov_data <- eig(avg_cov_data)
  cpc_error <- map_dbl(1:p, function(i){
      mean(map_dbl(reconstruct_cov_data, ~sum(.[i,-i]^2/(eig_avg_cov_data[i]*eig_avg_cov_data[-i])))) * ns / (p-1)
    })
  common_eigs <- map(reconstruct_cov_data, ~diag(.)[order(cpc_error)]) %>% unlist %>% matrix(p)
  gamma <- gamma[,order(cpc_error)]
  summary_mat <- data.frame(n_cpc = c(0:(p-2),p), 
                            error = sort(cpc_error),
                            quantile_l = NA,
                            quantile_u = NA)
  for(k in 0: (p-2)){
    # k is the number of cpc
    null_distribution <- map_dbl(1:m, function(j){
      if(k == 0){
        subject_eigs <- map(cov_data, ~eig(.)) %>% unlist %>% matrix(p)
        sim_cov_sqrt <- map(1:n, function(i){
          u <- svd(matrix(rnorm(p^2), p, p))$u
          d <- diag(sqrt(subject_eigs[,i]))
          u %*% d %*% t(u)
        })
        sim_data <- map(sim_cov_sqrt, ~ . %*% cov(matrix(rnorm(ns * p), ncol = p)) %*% .)
        avg_sim_data <- avg_mat(sim_data)
        pc_avg_sim_data <- svd(avg_sim_data)$u
        reconstruct_sim_data <- map(sim_data, ~t(pc_avg_sim_data) %*% . %*% pc_avg_sim_data)
        eig_avg_sim_data <- eig(avg_sim_data)
        map_dbl(1:p, function(i){mean(map_dbl(reconstruct_sim_data, 
                         ~sum(.[i,-i]^2/(eig_avg_sim_data[i]*eig_avg_sim_data[-i])))) * ns / (p-1)}) %>% min
      }else{
        gamma_true <- gamma[,1:k]
        subject_eigs <- map(1:n, ~eig(cov_data[[.]] - gamma[,1:k] %*% diag(common_eigs[1:k,.],nrow = k) %*% 
                                        t(gamma[,1:k]))) %>% unlist %>% matrix(p)
        sim_cov_sqrt <- tryCatch(map(1:n, function(i) {
          temp <- matrix(rnorm(p*(p-k)), p, p-k)
          u_i <- gramSchmidt(cbind(gamma_true, temp))$Q[,(k+1):p]
          gamma_true %*% diag(sqrt(common_eigs[1:k,i]), nrow = k) %*% t(gamma_true) +
            u_i %*% diag(sqrt(subject_eigs[1:(p-k),i]), nrow = p-k) %*% t(u_i)
        }), error = function(c){return(NA)})
        sim_data <- map(sim_cov_sqrt, ~ . %*% cov(matrix(rnorm(ns * p), ncol = p)) %*% .)
        if(min(is.na(sim_data))){return(NA)}
        avg_sim_data <- avg_mat(sim_data)
        pc_avg_sim_data <- svd(avg_sim_data)$u
        reconstruct_sim_data <- map(sim_data, ~t(pc_avg_sim_data) %*% . %*% pc_avg_sim_data)
        eig_avg_sim_data <- eig(avg_sim_data)
        sim_error <- map_dbl(1:p, function(i){
          mean(map_dbl(reconstruct_sim_data, 
                       ~sum(.[i,-i]^2/(eig_avg_sim_data[i]*eig_avg_sim_data[-i])))
          ) * ns / (p-1)
        })
        sim_error[order(sim_error)[k+1]]
      }
    })
    
    summary_mat[k+1, 3:4] <- quantile(null_distribution, c(0.05,1), na.rm = T)
    if(summary_mat$error[k+1] > summary_mat$quantile_l[k+1]){
      cpc_index <- order(cpc_error)[1:k]
      if(k == 0){cpc_index <- NULL}
      break
    }else if(k == p-2){
      cpc_index <- order(cpc_error)[1:(p-1)]
    }
  }
  return(list(cpc_index = cpc_index, summary_mat = summary_mat, cpc = gamma[,cpc_index]))
}


avg_mat <- function(list_of_mat){
  sum_of_mat <- 0
  for(i in 1:length(list_of_mat)){
    sum_of_mat <- sum_of_mat + list_of_mat[[i]]
  }
  sum_of_mat/length(list_of_mat)
}