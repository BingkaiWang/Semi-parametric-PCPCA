library(tidyverse)
library(GGally)
library(cowplot)
library(magick)
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

########### Visualization of Algorithm 2 --------
load("tfMRI_MOTOR_RL/Data.RData")
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
png(file = "data-analysis/figures/taskSM-semi.png")
heatmap(abs(cpc_taskSM), Colv = NA, Rowv = NA, scale = "none", revC = T)
dev.off()

############ Visualization of CPCs in brain view --------
# Generating Node file for BrainNet Viewer
for(i in c(1:30)){
  cpci <- cpc_taskSM[,i]
  roi_cpci <- filter(ROI.info.v2, module %in% c("SM_Hand", "SM_Mouth")) %>%
    dplyr::select(x, y, z) %>%
    mutate(color = ifelse(cpci > 0, 2, 1),
           size = ifelse(abs(cpci) > 0.15, abs(cpci), 0),
           label = '-') %>%
    filter(size != 0)
  write.table(roi_cpci,
              file = paste0("data-analysis/figures/semi-CPC/brainview/taskSM", i, "_brainview.node"),
              col.names = F, row.names = F)
}

# Then use BriainNet Viewer in MATLAB to generate brain maps.
# Results are stored in data-analysis/figures/semi-CPC/brainview/.

############# visualize of CPCs in time course -----------
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
ns <- nrow(ts.mc[[1]]) - 5
n <- length(ts)
ts.mc.SM <- map(ts.mc, ~(.[-c(1:5),SM.indi]))
ts.mc.SM.filtered <- map(ts.mc.SM, function(d){
  for(j in 1: p){d[,j] <- arma(d[,j], order = c(1,1))$residual}
  d[-1,]
})
ts.mc.SM.scaled <- map(ts.mc.SM, ~scale(.))

# get time brackets for tasks
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

for(j in 1:30){
  ggplot() + 
    geom_line(aes(x = (1:ns), y = abs(avg_mat(ts.mc.SM.scaled) %*% cpc_taskSM[,j])), size = 1) +
    geom_rect(aes(xmin = tongue_start-5, xmax = tongue_end-5, fill = "tongue"), ymin = -Inf, ymax = Inf, alpha = 0.2) +
    geom_rect(aes(xmin = lf_start-5, xmax = lf_end-5, fill = "left foot"), ymin = -Inf, ymax = Inf, alpha = 0.2) +
    geom_rect(aes(xmin = rf_start-5, xmax = rf_end-5, fill = "right foot"), ymin = -Inf, ymax = Inf, alpha = 0.2) +
    geom_rect(aes(xmin = lh_start-5, xmax = lh_end-5, fill = "left hand"), ymin = -Inf, ymax = Inf, alpha = 0.2) +
    geom_rect(aes(xmin = rh_start-5, xmax = rh_end-5, fill = "right hand"), ymin = -Inf, ymax = Inf, alpha = 0.2) +
    xlab("scan") + ylab("") + 
    theme(text = element_text(size=24),
          legend.text = element_text(size=20), 
          legend.position = "bottom",
          legend.spacing.x = unit(0.5, 'cm')
          ) +
    guides(fill = guide_legend(title = "Tasks", title.position = "left"))
   ggsave(paste0("data-analysis/figures/semi-CPC/timecourse/taskSM-tc-cpc", j, ".jpg"))
}


######### combining brainview and timecourse -----------------
taskSM_cor <- map(ts.mc.SM.filtered, ~cor(.))
avg_taskSM_cor <- avg_mat(taskSM_cor)
pc_taskSM_cor <- svd(avg_taskSM_cor)$u
eig_avg_taskSM_cor <- eig(avg_taskSM_cor)
reconstruct_taskSM_cor <- map(taskSM_cor, ~t(pc_taskSM_cor) %*% . %*% pc_taskSM_cor)
cpc_error <- map_dbl(1:p, function(i){
  mean(map_dbl(reconstruct_taskSM_cor,~sum(.[i,-i]^2/(eig_avg_taskSM_cor[i]*eig_avg_taskSM_cor[-i])))) * ns / (p-1)
})
var_explained <- eig_avg_taskSM_cor/sum(eig_avg_taskSM_cor)
var_explained <- paste0(round(var_explained[order(cpc_error)], 3) * 100, '%')
for(j in 1:30){
  path1 <- paste0("data-analysis/figures/semi-CPC/timecourse/taskSM-tc-cpc", j, ".jpg")
  path2 <- paste0("data-analysis/figures/semi-CPC/brainview/taskSM-cpc", j, ".jpg")
  title <- ggdraw() + draw_label(paste0("CPC ", j, ", variance explained: ", var_explained[j]), fontface = 'bold', size = 20)
  plot_tmp <- plot_grid(ggdraw() + draw_image(path1),
                        ggdraw() + draw_image(path2, scale = 0.9),
                        labels = NULL, ncol = 2)
  # plot_tmp <- plot_grid(title, plot_tmp, ncol = 1, rel_heights = c(0.1, 1))
  save_plot(paste0("data-analysis/figures/semi-CPC/taskSM-cpc", j, ".jpg"), plot_tmp,
            base_height = 5, base_aspect_ratio = 2.5)
}

########## visualization of nCPC ----------------------
ncpc_taskSM <- pc_taskSM_cor[,which(apply(abs(t(cpc_taskSM) %*% pc_taskSM_cor), 2, max) < 0.01)]
for(i in 1:5){
  cpci <- ncpc_taskSM[,i]
  roi_cpci <- filter(ROI.info.v2, module %in% c("SM_Hand", "SM_Mouth")) %>%
    dplyr::select(x, y, z) %>%
    mutate(color = ifelse(cpci > 0, 2, 1),
           size = ifelse(abs(cpci) > 0.15, abs(cpci), 0),
           label = '-') %>%
    filter(size != 0)
  write.table(roi_cpci,
              file = paste0("data-analysis/figures/semi-nCPC/brainview/ncpc", i, "_brainview.node"),
              col.names = F, row.names = F)
}
for(j in 1:ncol(ncpc_taskSM)){
  ggplot() + 
    geom_line(aes(x = (1:ns), y = abs(avg_mat(ts.mc.SM.scaled) %*% ncpc_taskSM[,j])), size = 1) +
    geom_rect(aes(xmin = tongue_start-5, xmax = tongue_end-5, fill = "tongue"), ymin = -Inf, ymax = Inf, alpha = 0.2) +
    geom_rect(aes(xmin = lf_start-5, xmax = lf_end-5, fill = "left foot"), ymin = -Inf, ymax = Inf, alpha = 0.2) +
    geom_rect(aes(xmin = rf_start-5, xmax = rf_end-5, fill = "right foot"), ymin = -Inf, ymax = Inf, alpha = 0.2) +
    geom_rect(aes(xmin = lh_start-5, xmax = lh_end-5, fill = "left hand"), ymin = -Inf, ymax = Inf, alpha = 0.2) +
    geom_rect(aes(xmin = rh_start-5, xmax = rh_end-5, fill = "right hand"), ymin = -Inf, ymax = Inf, alpha = 0.2) +
    xlab("scan") + ylab("") + 
    theme(text = element_text(size=24),
          legend.text = element_text(size=20), 
          legend.position = "bottom",
          legend.spacing.x = unit(0.5, 'cm')
    ) +
    guides(fill = guide_legend(title = "Tasks", title.position = "left"))
  ggsave(paste0("data-analysis/figures/semi-nCPC/timecourse/taskSM-tc-ncpc", j, ".jpg"))
}
for(j in 1:ncol(ncpc_taskSM)){
  path1 <- paste0("data-analysis/figures/semi-nCPC/timecourse/taskSM-tc-ncpc", j, ".jpg")
  path2 <- paste0("data-analysis/figures/semi-nCPC/brainview/ncpc", j, "_brainview.jpg")
  plot_tmp <- plot_grid(ggdraw() + draw_image(path1),
                        ggdraw() + draw_image(path2, scale = 0.9),
                        labels = NULL, ncol = 2)
  save_plot(paste0("data-analysis/figures/semi-nCPC/taskSM-cmb-ncpc", j, ".jpg"), plot_tmp,
            base_height = 5, base_aspect_ratio = 2.5)
}
