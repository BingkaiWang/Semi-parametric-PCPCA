set.seed(123)
library(tidyverse)
library(multigroup)
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

# correlation analysis on SM network using Flury's method
taskSM_stack <- NULL
taskSM_group <- NULL
for(j in 1:n){
  taskSM_stack <- rbind(taskSM_stack, ts.mc.SM.filtered[[j]])
  taskSM_group <- c(taskSM_group, rep(j, ns))
}
flury_result <- FCPCA(Data = taskSM_stack, Group = taskSM_group, graph = F)
pc_taskSM_flury <- flury_result$loadings.common
colnames(pc_taskSM_flury) <- paste("CPC", 1:35)
rownames(pc_taskSM_flury) <- paste("node", 1:35)
png(file = "data-analysis/figures/taskSM-flury.png")
heatmap(abs(pc_taskSM_flury), Colv = NA, Rowv = NA, scale = "none", revC = T)
dev.off()

