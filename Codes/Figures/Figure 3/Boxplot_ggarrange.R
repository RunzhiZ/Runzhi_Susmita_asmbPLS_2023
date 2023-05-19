library(ggplot2)
library(dplyr)
library(ggpubr)

dist <- "lognormal"
q = 200
p = 200

path <- paste0("C:/OneDrive/OneDrive - University of Florida/My research/asmbPLS/PLS/Simulation Hiper asmbPLS/q.", q, ".p.", p, "/", dist, "_mean_imputation/Boxplots")
setwd(path)

load(paste0("Results_for_plot_q.", q, ".p.", p, ".RData"))

boxplot_generate <- function(data, c, e) {
  boxplot <- data %>% 
    filter(Censoring_rate == c & Effect_scale == e) %>%
    mutate(Methods = ifelse((Method == "asmbPLS_1")|(Method == "asmbPLS_cv"), "asmbPLS", "Others")) %>%
    ggplot(aes(x = Method, y = MSE_scaled, fill = Methods)) + 
    #geom_jitter(width=0.1) +
    scale_fill_manual(values=c("grey", "white")) +
    theme(legend.position = "none") + guides(fill = "none") +
    geom_boxplot() + xlab("") + ylab("MSE (scaled)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 60, hjust=1)) +
    facet_grid(Noise_to_signal_ratio~Beta) + 
    theme(strip.text.x = element_text(size = 11, color = "black"),
          strip.text.y = element_text(size = 11, color = "black"),
          axis.text.x = element_text(size = 10, color = "black"),
          axis.text.y = element_text(size = 10, color = "black"),
          axis.title.x = element_text(size = 12, color = "black"),
          axis.title.y = element_text(size = 12, color = "black"))
  return(boxplot)
}

#### 0.5 ####
c = 0.1
e = 0.5
boxplot_0.5_0.1 <- boxplot_generate(Results_all_for_plot, c, e)

c = 0.3
e = 0.5
boxplot_0.5_0.3 <- boxplot_generate(Results_all_for_plot, c, e)

c = 0.5
e = 0.5
boxplot_0.5_0.5 <- boxplot_generate(Results_all_for_plot, c, e)

c = 0.7
e = 0.5
boxplot_0.5_0.7 <- boxplot_generate(Results_all_for_plot, c, e)

boxplot_0.5 <- ggarrange(boxplot_0.5_0.1, boxplot_0.5_0.3, boxplot_0.5_0.5, boxplot_0.5_0.7, 
                          nrow = 2, ncol = 2,  labels = c("A", "B", "C", "D"), common.legend = T)
jpeg(file = paste0("boxplot_0.5.jpeg"), width = 3000, height = 2500, res = 210, quality = 1200)
boxplot_0.5
dev.off()

#### 2 ####
c = 0.1
e = 2
boxplot_2_0.1 <- boxplot_generate(Results_all_for_plot, c, e)

c = 0.3
e = 2
boxplot_2_0.3 <- boxplot_generate(Results_all_for_plot, c, e)

c = 0.5
e = 2
boxplot_2_0.5 <- boxplot_generate(Results_all_for_plot, c, e)

c = 0.7
e = 2
boxplot_2_0.7 <- boxplot_generate(Results_all_for_plot, c, e)

boxplot_2 <- ggarrange(boxplot_2_0.1, boxplot_2_0.3, boxplot_2_0.5, boxplot_2_0.7, 
                         nrow = 2, ncol = 2,  labels = c("A", "B", "C", "D"), common.legend = T)
jpeg(file = paste0("boxplot_2.jpeg"), width = 3000, height = 2500, res = 210, quality = 1200)
boxplot_2
dev.off()

########################################## asmbPLS ###########################################
load("Results_for_plot_asmbPLS.RData")
c = 0.1
e = 2
boxplot_asmbPLS_2_0.1 <- boxplot_generate(Results_all_for_plot, c, e)

c = 0.3
e = 2
boxplot_asmbPLS_2_0.3 <- boxplot_generate(Results_all_for_plot, c, e)

c = 0.5
e = 2
boxplot_asmbPLS_2_0.5 <- boxplot_generate(Results_all_for_plot, c, e)

c = 0.7
e = 2
boxplot_asmbPLS_2_0.7 <- boxplot_generate(Results_all_for_plot, c, e)
boxplot_asmbPLS_2 <- ggarrange(boxplot_asmbPLS_2_0.1, boxplot_asmbPLS_2_0.3, boxplot_asmbPLS_2_0.5, boxplot_asmbPLS_2_0.7, 
                       nrow = 2, ncol = 2,  labels = c("A", "B", "C", "D"), common.legend = T)
jpeg(file = paste0("boxplot_asmbPLS_2.jpeg"), width = 3000, height = 2500, res = 210, quality = 1200)
boxplot_asmbPLS_2
dev.off()
