# environment setup ----------
set.seed(123)
library(tidyverse)
library(nlme)
library(tseries)
library(pracma) # eigenvalue of a matrix; gramschimdt
library(MASS)
library(Matrix)
avg_mat <- function(list_of_mat){
  sum_of_mat <- 0
  for(i in 1:length(list_of_mat)){
    sum_of_mat <- sum_of_mat + list_of_mat[[i]]
  }
  sum_of_mat/length(list_of_mat)
}

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
SM.indi <- which(ROI.info.v2$module %in% c("SM_Hand", "SM_Mouth"))
p <- length(SM.indi)
ns <- nrow(ts.mc[[1]]) - 6
n <- length(ts)
ts.mc.SM <- map(ts.mc, ~(.[-(1:5),SM.indi]))
ts.mc.SM.filtered <- map(ts.mc.SM, function(d){
  for(j in 1: p){d[,j] <- arma(d[,j], order = c(1,1))$residual}
  d[-1,]
})

# correlation analysis on SM network -----------
k <- ncol(ts.mc.SM.filtered[[1]])
taskSM_cor <- map(ts.mc.SM.filtered, ~cor(.))
subject_pc <- map(taskSM_cor, ~(svd(.)$v)[,1:k])
subject_pc_stack <- NULL
for(j in 1:n){
  subject_pc_stack <- cbind(subject_pc_stack, subject_pc[[j]])
}
pc_taskSM_pvd <- svd(subject_pc_stack)$u[,1:k]
png(file = "data-analysis/figures/taskSM-pvd.png")
heatmap(abs(pc_taskSM_pvd), Colv = NA, Rowv = NA, scale = "none", revC = T)
dev.off()


# visualization
ts.mc.SM.scaled <- map(ts.mc.SM, ~scale(.))
tongue_start <- map(onset, ~(.[,4])) %>% unlist %>% matrix(20) %>% .[19:20, ] %>% apply(1, mean)
tongue_duration <- map(onset, ~(.[,5])) %>% unlist %>% matrix(20) %>% .[19:20, ] %>% apply(1, mean)
tongue_end <- tongue_start + tongue_duration
lf_start <- map(onset, ~(.[,4])) %>% unlist %>% matrix(20) %>% .[11:12, ] %>% apply(1, mean)
lf_duration <- map(onset, ~(.[,5])) %>% unlist %>% matrix(20) %>% .[11:12, ] %>% apply(1, mean)
lf_end <- lf_start + lf_duration
rf_start <- map(onset, ~(.[,4])) %>% unlist %>% matrix(20) %>% .[13:14, ] %>% apply(1, mean)
rf_duration <- map(onset, ~(.[,5])) %>% unlist %>% matrix(20) %>% .[13:14, ] %>% apply(1, mean)
rf_end <- rf_start + rf_duration
lh_start <- map(onset, ~(.[,4])) %>% unlist %>% matrix(20) %>% .[15:16, ] %>% apply(1, mean)
lh_duration <- map(onset, ~(.[,5])) %>% unlist %>% matrix(20) %>% .[15:16, ] %>% apply(1, mean)
lh_end <- lh_start + lh_duration
rh_start <- map(onset, ~(.[,4])) %>% unlist %>% matrix(20) %>% .[17:18, ] %>% apply(1, mean)
rh_duration <- map(onset, ~(.[,5])) %>% unlist %>% matrix(20) %>% .[17:18, ] %>% apply(1, mean)
rh_end <- rh_start + rh_duration
for(j in 1:35){
  ggplot() + 
    geom_line(aes(x = (1:(ns+1)), y = abs(avg_mat(ts.mc.SM.scaled) %*% pc_taskSM_pvd[,j])), size = 1) +
    geom_rect(aes(xmin = tongue_start-5, xmax = tongue_end-5, fill = "tongue"), ymin = -Inf, ymax = Inf, alpha = 0.2) +
    geom_rect(aes(xmin = lf_start-5, xmax = lf_end-5, fill = "left foot"), ymin = -Inf, ymax = Inf, alpha = 0.2) +
    geom_rect(aes(xmin = rf_start-5, xmax = rf_end-5, fill = "right foot"), ymin = -Inf, ymax = Inf, alpha = 0.2) +
    geom_rect(aes(xmin = lh_start-5, xmax = lh_end-5, fill = "left hand"), ymin = -Inf, ymax = Inf, alpha = 0.2) +
    geom_rect(aes(xmin = rh_start-5, xmax = rh_end-5, fill = "right hand"), ymin = -Inf, ymax = Inf, alpha = 0.2) +
    xlab("scan") + ylab("") +
    theme(text = element_text(size=24),
          legend.text = element_text(size=24), 
          legend.position = "bottom",
          legend.spacing.x = unit(0.5, 'cm')
    ) +
    guides(fill = guide_legend(title = "Tasks", title.position = "left"))
  ggsave(paste0("data-analysis/figures/PVD/taskSM-tc-pvd-cpc", j, ".jpg"), width = 12, height = 7)
}
