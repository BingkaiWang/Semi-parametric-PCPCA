# environment setup ----------
set.seed(1234)
library(multigroup)
library(nlme)
library(tseries)
library(tidyverse)
source("simulations/pcpc_test.R") # loading required functions and packages 
# parallelism
library(foreach)
library(doParallel)
library(tseries)
cl <- makeCluster(16)
registerDoParallel(cl)
m <- 1000

# data import and preprocessing ------------
load("tfMRI_MOTOR_RL/Data.RData")
# motion correction and temporal decorrelation
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
ns <- nrow(ts.mc[[1]]) - 6
n <- length(ts)
p <- ncol(ts.mc[[1]])
ts.mc.filtered <- map(ts.mc, function(d){
  d <- d[-c(1:5),]
  for(j in 1: p){d[,j] <- arma(d[,j], order = c(1,1))$residual}
  d[-1,]
})

# auto-correlation comparison
lag.max <- 25
ac.before <- map(ts.mc, function(d){
  d <- d[-c(1:5),]
  map(1:ncol(d), ~acf(d[,.], lag.max = lag.max, plot = F)$acf) %>% unlist %>% matrix(nrow = lag.max+1)
}) %>% do.call("cbind", .) %>% .[-1,]
ac.after <- map(ts.mc.filtered, function(d){
  d <- d[-c(1:5),]
  map(1:ncol(d), ~acf(d[,.], lag.max = lag.max, plot = F)$acf) %>% unlist %>% matrix(nrow = lag.max+1)
}) %>% do.call("cbind", .) %>% .[-1,]

plot.ac <- data.frame(Lag = factor(rep(1:lag.max, 2)), y = c(rowMeans(ac.before), rowMeans(ac.after)),
                      sd = c(apply(ac.before,1,sd), apply(ac.after,1,sd)),
                      y_lower = c(map_dbl(1:lag.max, ~quantile(ac.before[.,], 0.05)), map_dbl(1:lag.max, ~quantile(ac.after[.,], 0.05))),
                      y_upper = c(map_dbl(1:lag.max, ~quantile(ac.before[.,], 0.95)), map_dbl(1:lag.max, ~quantile(ac.after[.,], 0.95))),
                      Label = rep(c("Before", "After"), each = lag.max)
                      )
ggplot(filter(plot.ac, Lag %in% 1:15)) +
  geom_point(aes(x=Lag, y = y, color = Label), size = 4, position = position_dodge(width = 0.7)) +
  geom_errorbar(aes(x = Lag,  ymin = y_lower, ymax = y_upper, color = Label), size = 1.2, width = 0.7, position = position_dodge(width = 0.7)) +
  theme_bw() +
  theme(text = element_text(size = 35), 
        legend.position = c(0.85,0.85), 
        legend.text = element_text(size = 30), 
        legend.key.size = unit(2, "cm"), 
        legend.title = element_blank()) + 
  scale_y_continuous(name="Autocorrelation", limits=c(-0.25, 0.75))
ggsave("auto-correlation.png", width = 18, height =10)

# correlation analysis on all regions -----------
task_cor <- map(ts.mc.filtered, ~cor(.))
avg_task_cor <- avg_mat(task_cor)
pc_task_cor <- svd(avg_task_cor)$u

# find ordering of CPC ------------------
reconstruct_task_cor <- map(task_cor, ~t(pc_task_cor) %*% . %*% pc_task_cor)
eig_avg_task_cor <- eig(avg_task_cor)
cpc_error <- map_dbl(1:p, function(i){
  mean(map_dbl(reconstruct_task_cor, ~sum(.[i,-i]^2/(eig_avg_task_cor[i]*eig_avg_task_cor[-i])))) * ns / (p-1)
})
common_eigs <- map(reconstruct_task_cor, ~diag(.)[order(cpc_error)]) %>% 
  unlist %>% matrix(p)
pc_task_cor <- pc_task_cor[,order(cpc_error)]
summary_mat <- data.frame(n_cpc = c(0:(p-2),p), 
                          error = sort(cpc_error),
                          quantile_l = NA,
                          quantile_u = NA)
packages = c('tidyverse','pracma', 'MASS', 'Matrix')
for(k in 0: (p-2)){
  # k is the number of cpc
  null_distribution <- foreach(j=1:m, .combine = cbind, .packages = packages) %dopar% {
    if(k == 0){
      subject_eigs <- map(task_cor, ~eig(.)) %>% unlist %>% matrix(p)
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
      map_dbl(1:p, function(i){
        mean(map_dbl(reconstruct_sim_data,~sum(.[i,-i]^2/(eig_avg_sim_data[i]*eig_avg_sim_data[-i])))) * ns / (p-1)
      }) %>% min
    }else{
      gamma_true <- pc_task_cor[,1:k]
      subject_eigs <- map(1:n, ~eig(task_cor[[.]] - gamma_true %*% diag(common_eigs[1:k,.],nrow = k) %*% 
                                      t(gamma_true))) %>% unlist %>% matrix(p)
      sim_cov_sqrt <- map(1:n, function(i) {
        temp <- matrix(rnorm(p*(p-k)), p, p-k)
        while (rankMatrix(cbind(gamma_true, temp)) != p) {
          temp <- matrix(rnorm(p*(p-k)), p, p-k)
        }
        u_i <- gramSchmidt(cbind(gamma_true, temp))$Q[,(k+1):p]
        gamma_true %*% diag(sqrt(common_eigs[1:k,i]), nrow = k) %*% t(gamma_true) +
          u_i %*% diag(sqrt(subject_eigs[1:(p-k),i]), nrow = p-k) %*% t(u_i)
      })
      sim_data <- map(sim_cov_sqrt, ~ . %*% cov(matrix(rnorm(ns * p), ncol = p)) %*% .)
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
  write.table(null_distribution, file = paste0("testing-k=",k,".txt"), sep="\t", row.names=F)
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

saveRDS(list(cpc_index = cpc_index, summary_mat = summary_mat, cpc = pc_task_cor[,cpc_index]), "step2_task_RL.rds")



# visualization -------------
# Normality test of data
library(MVN)
d <- NULL
for(i in length(ts.mc.filtered)){
  d <- rbind(d, ts.mc.filtered[[i]])
}
colnames(d) <- paste0("V",1:264)
d <- scale(d, center = T, scale = F)
mvn(data = d,mvnTest = "mardia")$multivariateNormality
mvn(data = d,mvnTest = "hz")$multivariateNormality
mvn(data = d,mvnTest = "royston")$multivariateNormality
mvn(data = d,mvnTest = "dh")$multivariateNormality
mvn(data = d,mvnTest = "energy")$multivariateNormality
# 
# 
task_cor <- map(ts.mc.filtered, ~cor(.))
avg_task_cor <- avg_mat(task_cor)
pc_task_cor <- svd(avg_task_cor)$u
reconstruct_task_cor <- map(task_cor, ~t(pc_task_cor) %*% . %*% pc_task_cor)
eig_avg_task_cor <- eig(avg_task_cor)
cpc_error <- map_dbl(1:p, function(i){
  mean(map_dbl(reconstruct_task_cor, ~sum(.[i,-i]^2/(eig_avg_task_cor[i]*eig_avg_task_cor[-i])))) * ns / (p-1)
})
cpc_hat <- pc_task_cor[,order(cpc_error)]

CPC_associated_networks <- list(1:264)
for(j in 1:264){
  ne <- ROI.info.v2$module[abs(cpc_hat[,j]) > 0.1] %>% table
  CPC_associated_networks[[j]] <- names(ne)[ne/sum(ne) > 0.25]
}
CPC_associated_networks[1:190] %>% unlist %>% table

for(j in 1:264){
  ggplot(data.frame(x = cpc_hat[,j], label = ROI.info.v2$module)) + 
    geom_point(aes(x=1:264, y = abs(x * (abs(x) > 0.1)), shape = label, color = label), size = 3) +
    scale_shape_manual(values=1:nlevels(ROI.info.v2$module))
  ggsave(paste0("all-region-loading/cpc", j, ".jpg"))
}
SM.indi <- which(ROI.info.v2$module %in% c("SM_Hand", "SM_Mouth"))

xx <-  pc_task_cor[SM.indi,order(cpc_error)]
xx <- apply(xx, 2, function(x){x/sqrt(sum(x^2))})
apply(abs(t(cpc_taskSM) %*% xx),1, max)
