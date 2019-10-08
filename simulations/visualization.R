# simulation 1 visualziation
library(tidyverse)
library(cowplot)
sim1 <- readRDS("simulations/simulation1.rds")
sim1_n <- t(sim1$sim1_n)
sim1_ns <- t(sim1$sim1_ns)

sim1_n_plot <- data.frame(err =  as.vector(sim1_n),
                          n = rep(rep(c(50, 100,  500, 1000), each = 1000),2),
                          label = rep(c("CPC estimate", "non-CPC estimate"), each = 4000)
                          )

sim1_n_figure <- ggplot(sim1_n_plot) +
  geom_violin(aes(x = factor(n), y = err, color = label), position = position_dodge(width = 0.3)) +
  geom_hline(aes(yintercept = 1/50), alpha = 0.5) +
  geom_hline(aes(yintercept = 0.3027), alpha = 0.5, linetype = "dashed") +
  xlab("n") + ylab("Deviation from Commonality Metric") + 
  theme(legend.position = c(0,0.95),
        legend.title = element_blank(),
        text = element_text(size = 12)
  )


sim1_ns_plot <- data.frame(err =  as.vector(sim1_ns),
                          ns = rep(rep(c(50, 100,  500, 1000), each = 1000),2),
                          label = rep(c("CPC estimate", "non-CPC estimate"), each = 4000)
)

sim1_ns_figure <- ggplot(sim1_ns_plot) +
  geom_violin(aes(x = factor(ns), y = err, color = label), position = position_dodge(width = 0.3)) +
  geom_hline(aes(yintercept = 0), alpha = 0.5) +
  geom_hline(aes(yintercept = 0.1073), alpha = 0.5, linetype = "dashed") +
  xlab("T") + ylab("Deviation from Commonality Metric") + 
  theme(legend.position = c(0,0.95),
        legend.title = element_blank(),
        text = element_text(size = 12)
  )

sim1_plot <- plot_grid(sim1_n_figure, sim1_ns_figure, ncol = 2)
save_plot("simulation1.png", sim1_plot, base_aspect_ratio = 2)

# simulation 3 visualization
sim3.1 <- readRDS("simulations/simulation3-1.rds")
sim3.1 <- sim3.1[4:6,]
apply(sim3.1, 1, mean)
sim3.2 <- readRDS("simulations/simulation3-2.rds")
sim3.2 <- sim3.2[4:6,]
apply(sim3.2, 1, mean)