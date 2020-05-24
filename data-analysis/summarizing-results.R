library(tidyverse)
library(GGally)

# SM visualization -----
step2_taskSM <- readRDS("data-analysis/step2_taskSM_RL.rds")
cpc_taskSM <- step2_taskSM$cpc
plot_data_taskSM <- step2_taskSM$summary_mat

ggplot(plot_data_taskSM) +
  geom_point(aes(x = factor(n_cpc), y = error, shape = "Sensorimotor")) +
  geom_errorbar(aes(x = factor(n_cpc), ymin = quantile_l, ymax = quantile_u, linetype = "Reference C.I."), width = 0.5) +
  theme_bw() +
  theme(text = element_text(size = 12),
        legend.title=element_blank(),
        legend.justification=c(1,1), legend.position=c(1,1),
        legend.background=element_blank(),
        legend.text = element_text(size = 12)) +
  xlab("Number of CPCs to test") + ylab("Error Measure")
ggsave("data-analysis/figures/taskSM-step2.png")

# region analysis -----
load("tfMRI_MOTOR_RL/Data.RData")
source("simulations/pcpc_test.R")
thresholding <- function(x, thresh){
  ifelse(abs(x) > thresh, x, 0)
}
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
# correlation analysis on SM network
SM.indi <- which(ROI.info.v2$module %in% c("SM_Hand", "SM_Mouth"))
p <- length(SM.indi)
ns <- nrow(ts.mc[[1]])
ts.mc.SM <- map(ts.mc, ~scale(.[-(1:3),SM.indi]))
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

ggplot() + 
  geom_line(aes(x = 1:(ns-3), y = abs(avg_mat(ts.mc.SM) %*% cpc_taskSM[,29]))) +
  geom_rect(aes(xmin = tongue_start, xmax = tongue_end, fill = "tongue"), ymin = -Inf, ymax = Inf, alpha = 0.2) +
  geom_rect(aes(xmin = lf_start, xmax = lf_end, fill = "left foot"), ymin = -Inf, ymax = Inf, alpha = 0.2) +
  geom_rect(aes(xmin = rf_start, xmax = rf_end, fill = "right foot"), ymin = -Inf, ymax = Inf, alpha = 0.2) +
  geom_rect(aes(xmin = lh_start, xmax = lh_end, fill = "left hand"), ymin = -Inf, ymax = Inf, alpha = 0.2) +
  geom_rect(aes(xmin = rh_start, xmax = rh_end, fill = "right hand"), ymin = -Inf, ymax = Inf, alpha = 0.2) +
  xlab("scan") + ylab("") +
  theme(legend.text = element_text(size=16), legend.title = element_blank())
ggsave(paste0("Motor-task/motor-task-plot-RL/cpc_timecourse/taskSM-tc-cpc", j, ".jpg"))
