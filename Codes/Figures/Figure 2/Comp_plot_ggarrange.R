library(ggplot2)
library(dplyr)
library(ggpubr)
comp_summary_lognormal <- NULL
comp_summary_weibull <- NULL

dist <- "lognormal"
q = 200
p = 200

for (c in c("00", "01", "03", "05", "06", "07")) {
  for (r in c("00", "01", "02", "05", "10")) {
    for (effect in c("05", "20")) {
      setwd(paste0("C:/OneDrive/OneDrive - University of Florida/My research/asmbPLS/PLS/Simulation Hiper asmbPLS/q.", q, ".p.", p, "/", dist, "_mean_imputation/Results_V4_plus"))
      eval(parse(text = paste0("load(\"Results_c", c, "_r", r, "_e", effect, "_", dist, ".RData\")")))
      for (n in 1:6) {
        temp <- as.data.frame(matrix(NA, nrow = 5, ncol = 7))
        colnames(temp) <- c("Distribution", "Censoring_rate", "Noise_to_signal_ratio", "Effect_scale",
                            "Beta", "Component", "Value")
        c_input <- c("0", "0.1", "0.3", "0.5", "0.6", "0.7")[which(c("00", "01", "03", "05", "06", "07") %in% c)]
        r_input <- c("0", "0.1", "0.2", "0.5", "1")[which(c("00", "01", "02", "05", "10") %in% r)]
        effect_input <- c("0.5", "2")[which(c("05", "20") %in% effect)]
        temp$Censoring_rate <- c_input
        temp$Noise_to_signal_ratio <- r_input
        temp$Effect_scale <- effect_input
        temp$Beta <- paste0("beta ", n)
        temp$Component <- c("comp_1", "comp_2", "comp_3", "comp_4", "comp_5")
        asmbPLS_results <- results_asmbPLS_0.7_0.7_predict_list[[n]]
        comp_temp <- apply(asmbPLS_results, 2, which.min)
        temp$Value <- c(length(which(comp_temp == 1)),
                        length(which(comp_temp == 2)),
                        length(which(comp_temp == 3)),
                        length(which(comp_temp == 4)),
                        length(which(comp_temp == 5)))
        comp_summary_lognormal <- rbind(comp_summary_lognormal, temp)
      }
    }
  }
}

dist <- "weibull"

for (c in c("00", "01", "03", "05", "06", "07")) {
  for (r in c("00", "01", "02", "05", "10")) {
    for (effect in c("05", "20")) {
      setwd(paste0("C:/OneDrive/OneDrive - University of Florida/My research/asmbPLS/PLS/Simulation Hiper asmbPLS/q.", q, ".p.", p, "/", dist, "_mean_imputation/Results_V4"))
      eval(parse(text = paste0("load(\"Results_c", c, "_r", r, "_e", effect, "_", dist, ".RData\")")))
      for (n in 1:6) {
        temp <- as.data.frame(matrix(NA, nrow = 5, ncol = 7))
        colnames(temp) <- c("Distribution", "Censoring_rate", "Noise_to_signal_ratio", "Effect_scale",
                            "Beta", "Component", "Value")
        c_input <- c("0", "0.1", "0.3", "0.5", "0.6", "0.7")[which(c("00", "01", "03", "05", "06", "07") %in% c)]
        r_input <- c("0", "0.1", "0.2", "0.5", "1")[which(c("00", "01", "02", "05", "10") %in% r)]
        effect_input <- c("0.5", "2")[which(c("05", "20") %in% effect)]
        temp$Censoring_rate <- c_input
        temp$Noise_to_signal_ratio <- r_input
        temp$Effect_scale <- effect_input
        temp$Beta <- paste0("beta ", n)
        temp$Component <- c("comp_1", "comp_2", "comp_3", "comp_4", "comp_5")
        asmbPLS_results <- results_asmbPLS_0.7_predict_list[[n]]
        comp_temp <- apply(asmbPLS_results, 2, which.min)
        temp$Value <- c(length(which(comp_temp == 1)),
                        length(which(comp_temp == 2)),
                        length(which(comp_temp == 3)),
                        length(which(comp_temp == 4)),
                        length(which(comp_temp == 5)))
        comp_summary_weibull <- rbind(comp_summary_weibull, temp)
      }
    }
  }
}

comp_plot_lognormal_0.5 <- comp_summary_lognormal %>% 
  filter(Effect_scale == 0.5) %>%
  filter((Censoring_rate == 0.1)|(Censoring_rate == 0.3)|(Censoring_rate == 0.5)|(Censoring_rate == 0.7)) %>%
  ggplot(aes(x = Beta, y = Value, fill = Component)) + 
  #geom_jitter(width=0.1) +
  theme(legend.position = "none") + 
  geom_bar(stat = 'identity') + xlab("") + ylab("Proportion") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  facet_grid(Noise_to_signal_ratio ~ Censoring_rate) + 
  theme(strip.text.x = element_text(size = 11, color = "black"),
        strip.text.y = element_text(size = 11, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"))

comp_plot_lognormal_2 <- comp_summary_lognormal %>% 
  filter(Effect_scale == 2) %>%
  filter((Censoring_rate == 0.1)|(Censoring_rate == 0.3)|(Censoring_rate == 0.5)|(Censoring_rate == 0.7)) %>%
  ggplot(aes(x = Beta, y = Value, fill = Component)) + 
  #geom_jitter(width=0.1) +
  theme(legend.position = "none") + 
  geom_bar(stat = 'identity') + xlab("") + ylab("Proportion") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  facet_grid(Noise_to_signal_ratio ~ Censoring_rate) + 
  theme(strip.text.x = element_text(size = 11, color = "black"),
        strip.text.y = element_text(size = 11, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"))


comp_plot_weibull_0.5 <- comp_summary_weibull %>% 
  filter(Effect_scale == 0.5) %>%
  filter((Censoring_rate == 0.1)|(Censoring_rate == 0.3)|(Censoring_rate == 0.5)|(Censoring_rate == 0.7)) %>%
  ggplot(aes(x = Beta, y = Value, fill = Component)) + 
  #geom_jitter(width=0.1) +
  theme(legend.position = "none") + 
  geom_bar(stat = 'identity') + xlab("") + ylab("Proportion") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  facet_grid(Noise_to_signal_ratio ~ Censoring_rate) + 
  theme(strip.text.x = element_text(size = 11, color = "black"),
        strip.text.y = element_text(size = 11, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"))

comp_plot_weibull_2 <- comp_summary_weibull %>% 
  filter(Effect_scale == 2) %>%
  filter((Censoring_rate == 0.1)|(Censoring_rate == 0.3)|(Censoring_rate == 0.5)|(Censoring_rate == 0.7)) %>%
  ggplot(aes(x = Beta, y = Value, fill = Component)) + 
  #geom_jitter(width=0.1) +
  theme(legend.position = "none") + 
  geom_bar(stat = 'identity') + xlab("") + ylab("Proportion") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  facet_grid(Noise_to_signal_ratio ~ Censoring_rate) + 
  theme(strip.text.x = element_text(size = 11, color = "black"),
        strip.text.y = element_text(size = 11, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"))

setwd(paste0("C:/OneDrive/OneDrive - University of Florida/My research/asmbPLS/PLS/Simulation Hiper asmbPLS/q.", q, ".p.", p))

comp_plot <- ggarrange(comp_plot_lognormal_0.5, comp_plot_lognormal_2, comp_plot_weibull_0.5, comp_plot_weibull_2, 
                       nrow = 2, ncol = 2,  labels = c("A", "B", "C", "D"), common.legend = T)

jpeg(file = "Comp_plot_ggarrange.jpeg", width = 3000, height = 2500, res = 250, quality = 1200)
comp_plot
dev.off()


