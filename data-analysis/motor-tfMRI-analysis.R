# environment setup ----------
set.seed(123)
library(multigroup)
library(nlme)
source("~/graphical_model/PCPCA/simulations/pcpc_test.R") # loading required functions and packages 
# parallelism
library(foreach)
library(doParallel)
cl <- makeCluster(32)
registerDoParallel(cl)
m <- 1000

# data import and preprocessing ------------
load("~/graphical_model/PCPCA/tfMRI_MOTOR_RL/Data.RData")
# motion correction
ts <- ts[ts.indi == 1]
motion <- motion[ts.indi == 1]
ts.mc <- ts
n <- length(ts)
for(i in 1:n){
  for(j in 1:ncol(ts[[i]])){
    fit.tmp <- lm(ts[[i]][,j]~motion[[i]])
    ts.mc[[i]][,j] <- fit.tmp$residuals
  }
}
SM.indi <- which(ROI.info.v2$module %in% c("SM_Hand", "SM_Mouth"))
p <- length(SM.indi)
ns <- nrow(ts.mc[[1]])
n <- length(ts)

# correlation analysis on SM network -----------
ts.mc.SM <- map(ts.mc, ~(.[-(1:5),SM.indi]))
taskSM_cor <- map(ts.mc.SM, ~cor(.))
avg_taskSM_cor <- avg_mat(taskSM_cor)
pc_taskSM_cor <- svd(avg_taskSM_cor)$u

# find ordering of CPC ------------------
reconstruct_taskSM_cor <- map(taskSM_cor, ~t(pc_taskSM_cor) %*% . %*% pc_taskSM_cor)
cpc_error <- map_dbl(1:p, function(i){
  map_dbl(reconstruct_taskSM_cor, ~pc_cor(., index = i)) %>% mean
})
eig_avg_taskSM_cor <- eig(avg_taskSM_cor)
cpc_error <- map_dbl(1:p, function(i){
  mean(map_dbl(reconstruct_taskSM_cor, ~sum(.[i,-i]^2/(eig_avg_taskSM_cor[i]*eig_avg_taskSM_cor[-i])))) * ns / (p-1)
})
reconstruct_eigenvalues <- map(reconstruct_taskSM_cor, ~diag(.)[order(cpc_error)]) %>% 
  unlist %>% matrix(p)
pc_taskSM_cor <- pc_taskSM_cor[,order(cpc_error)]
summary_mat <- data.frame(n_cpc = c(0:(p-2),p), 
                          error = sort(cpc_error),
                          quantile_l = NA,
                          quantile_u = NA)
packages = c('tidyverse','pracma', 'MASS', 'Matrix')
for(k in 0: (p-2)){
  # k is the number of cpc
  null_distribution <- foreach(j=1:m, .combine = cbind, .packages = packages) %dopar% {
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
      eig_avg_sim_data <- eig(avg_sim_data)
      map_dbl(1:p, function(i){
        mean(map_dbl(reconstruct_sim_data,~sum(.[i,-i]^2/(eig_avg_sim_data[i]*eig_avg_sim_data[-i])))) * ns / (p-1)
      }) %>% min
    }else{
      gamma_true <- pc_taskSM_cor[,1:k]
      sim_data <- map(1:n, function(i) {
        temp <- matrix(rnorm(p*(p-k)), p, p-k)
        while (rankMatrix(cbind(gamma_true, temp)) != p) {
          temp <- matrix(rnorm(p*(p-k)), p, p-k)
        }
        u_i <- gramSchmidt(cbind(gamma_true, temp))$Q[,(k+1):p]
        gamma_true %*% diag(reconstruct_eigenvalues[1:k,i], nrow = k) %*% t(gamma_true) +
          u_i %*% diag(reconstruct_eigenvalues[(k+1):p,i], nrow = p-k) %*% t(u_i)
      })
      sim_data <- map(sim_data, ~cov(mvrnorm(ns, rep(0,p), .)))
      avg_sim_data <- avg_mat(sim_data)
      pc_avg_sim_data <- svd(avg_sim_data)$u
      reconstruct_sim_data <- map(sim_data, ~t(pc_avg_sim_data) %*% . %*% pc_avg_sim_data)
      eig_avg_sim_data <- eig(avg_sim_data)
      sim_error <- map_dbl((1:p), function(i){
        mean(map_dbl(reconstruct_sim_data, ~sum(.[i,-i]^2/(eig_avg_sim_data[i]*eig_avg_sim_data[-i])))) * ns / (p-1)
      })
      sim_error[order(sim_error)[k+1]]
    }
  }
  null_distribution <- as.vector(null_distribution)
  summary_mat[k+1, 3:4] <- quantile(null_distribution, c(0.05,1))
  if(summary_mat$error[k+1] > summary_mat$quantile_l[k+1]){
    cpc_index <- order(cpc_error)[1:k]
    if(k == 0){cpc_index <- NULL}
    break
  }else if(k == p-2){
    cpc_index <- order(cpc_error)[1:(p-1)]
  }
}

saveRDS(list(cpc_index = cpc_index, summary_mat = summary_mat, cpc = pc_taskSM_cor[,cpc_index]), "step2_taskSM_RL.rds")

# normality test
library(MVN)
d <- NULL
for(i in length(ts.mc.SM)){
  d <- rbind(d, ts.mc.SM[[i]])
}
colnames(d) <- paste0("V",1:35)
d <- scale(d, center = T, scale = F)
mvn(data = d,mvnTest = "mardia")$multivariateNormality
mvn(data = d,mvnTest = "hz")$multivariateNormality
mvn(data = d,mvnTest = "royston")$multivariateNormality
mvn(data = d,mvnTest = "dh")$multivariateNormality
mvn(data = d,mvnTest = "energy")$multivariateNormality