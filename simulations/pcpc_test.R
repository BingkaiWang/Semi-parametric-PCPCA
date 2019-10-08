library(purrr)
library(dplyr)
library(pracma)
library(MASS)
library(Matrix)

pcpc_test <- function(gamma, cov_data, ns, m = 100, err_measure = 1){
  reconstruct_cov_data <- map(cov_data, ~t(gamma) %*% . %*% gamma)
  p <- ncol(gamma)
  n <- length(cov_data)
  avg_cov_data <- avg_mat(cov_data)
  eig_avg_cov_data <- eig(avg_cov_data)
  if(err_measure == 1){
    cpc_error <- map_dbl(1:p, function(i){
      mean(map_dbl(reconstruct_cov_data, 
                   ~sum(.[i,-i]^2/(eig_avg_cov_data[i]*eig_avg_cov_data[-i])))
      ) * ns / (p-1)
    })
  }else if(err_measure == 2){
    cpc_error <- map_dbl(1:p, function(i){
      map_dbl(reconstruct_cov_data, ~pc_cor(., index =i)) %>% mean
    })
  }else{
    stop()
  }
  reconstruct_eigenvalues <- map(reconstruct_cov_data, ~diag(.)[order(cpc_error)]) %>% 
    unlist %>% matrix(p)
  gamma <- gamma[,order(cpc_error)]
  summary_mat <- data.frame(n_cpc = c(0:(p-2),p), 
                            error = sort(cpc_error),
                            quantile_l = NA,
                            quantile_u = NA)
  for(k in 0: (p-2)){
    # k is the number of cpc
    null_distribution <- map_dbl(1:m, function(j){
      if(k == 0){
        sim_data <- map(1:n, function(i){
          u <- svd(matrix(rnorm(p^2), p, p))$u
          d <- diag(reconstruct_eigenvalues[,i])
          u %*% d %*% t(u)
        })
        sim_data <- map(sim_data, ~cov(mvrnorm(ns, rep(0,p), .)))
        avg_sim_data <- avg_mat(sim_data)
        pc_avg_sim_data <- svd(avg_sim_data)$u
        reconstruct_sim_data <- map(sim_data, ~t(pc_avg_sim_data) %*% . %*% pc_avg_sim_data)
        if(err_measure == 1){
          eig_avg_sim_data <- eig(avg_sim_data)
          map_dbl(1:p, function(i){
            mean(map_dbl(reconstruct_sim_data, 
                         ~sum(.[i,-i]^2/(eig_avg_sim_data[i]*eig_avg_sim_data[-i])))
            ) * ns / (p-1)
          }) %>% min
        }else{
          map_dbl(1:p, function(i){
            map_dbl(reconstruct_sim_data, ~pc_cor(., index =i)) %>% mean
          }) %>% min
        }
      }else{
        gamma_true <- gamma[,1:k]
        sim_data <- tryCatch(map(1:n, function(i) {
          temp <- matrix(rnorm(p*(p-k)), p, p-k)
          u_i <- gramSchmidt(cbind(gamma_true, temp))$Q[,(k+1):p]
          gamma_true %*% diag(reconstruct_eigenvalues[1:k,i], nrow = k) %*% t(gamma_true) +
            u_i %*% diag(reconstruct_eigenvalues[(k+1):p,i], nrow = p-k) %*% t(u_i)
        }), error = function(c){return(NA)})
        if(min(is.na(sim_data))){return(NA)}
        sim_data <- map(sim_data, ~cov(mvrnorm(ns, rep(0,p), .)))
        avg_sim_data <- avg_mat(sim_data)
        pc_avg_sim_data <- svd(avg_sim_data)$u
        reconstruct_sim_data <- map(sim_data, ~t(pc_avg_sim_data) %*% . %*% pc_avg_sim_data)
        if(err_measure == 1){
          eig_avg_sim_data <- eig(avg_sim_data)
          sim_error <- map_dbl(1:p, function(i){
            mean(map_dbl(reconstruct_sim_data, 
                         ~sum(.[i,-i]^2/(eig_avg_sim_data[i]*eig_avg_sim_data[-i])))
            ) * ns / (p-1)
          })
        } else{
          sim_error <- map_dbl(1:p, function(i){
            map_dbl(reconstruct_sim_data, ~pc_cor(., index =i)) %>% mean
          })
        }
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

pc_cor <- function(matrix, index){
  matrix <- cov2cor(matrix)
  if(length(index) == 1){
    sum((matrix[index,-index])^2)
  }else{
    diag(matrix) <- 0
    sum(matrix[index,]^2)
  }
}