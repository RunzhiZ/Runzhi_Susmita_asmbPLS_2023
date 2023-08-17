library(ggplot2)
library(dplyr)
dist <- "lognormal"
q = 200
p = 200
path <- paste0("C:/OneDrive/OneDrive - University of Florida/My research/asmbPLS/PLS/Simulation Hiper asmbPLS/q.", q, ".p.", p, "/", dist, "_mean_imputation/AUC_plots")
setwd(path)
load("Feature_selection_summary.RData")

Results_all_feature_selection_microbiome_summary$Method_effect <- paste0(Results_all_feature_selection_microbiome_summary$Method, "_",
                                                                         Results_all_feature_selection_microbiome_summary$Effect_scale)
Results_all_feature_selection_metabolome_summary$Method_effect <- paste0(Results_all_feature_selection_metabolome_summary$Method, "_",
                                                                         Results_all_feature_selection_metabolome_summary$Effect_scale)

#### microbiome ####
micro_sensitivity_plot <- Results_all_feature_selection_microbiome_summary %>% 
  filter(Censoring_rate %in% c("0.1", "0.3", "0.5", "0.7")) %>%
  ggplot(aes(x = Noise_to_signal_ratio, y = Sensitivity, color = Method, shape = Effect_scale, group = Method_effect)) + 
  theme(legend.position = "none") +
  geom_line() + geom_point() +
  xlab("Noise") + guides(fill = "none") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  facet_grid(Censoring_rate~Beta) + 
  theme(strip.text.x = element_text(size = 11, color = "black"),
        strip.text.y = element_text(size = 11, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 12, color = "black"))

micro_specificity_plot <- Results_all_feature_selection_microbiome_summary %>% 
  filter(Censoring_rate %in% c("0.1", "0.3", "0.5", "0.7")) %>%
  ggplot(aes(x = Noise_to_signal_ratio, y = Specificity, color = Method, shape = Effect_scale, group = Method_effect)) + 
  theme(legend.position = "none") +
  geom_line() + geom_point() +
  xlab("Noise") + guides(fill = "none") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  facet_grid(Censoring_rate~Beta) + 
  theme(strip.text.x = element_text(size = 11, color = "black"),
        strip.text.y = element_text(size = 11, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 12, color = "black"))

#### metabolome ####
meta_sensitivity_plot <- Results_all_feature_selection_metabolome_summary %>% 
  filter(Censoring_rate %in% c("0.1", "0.3", "0.5", "0.7")) %>%
  ggplot(aes(x = Noise_to_signal_ratio, y = Sensitivity, color = Method, shape = Effect_scale, group = Method_effect)) + 
  theme(legend.position = "none") +
  geom_line() + geom_point() +
  xlab("Noise") + guides(fill = "none") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  facet_grid(Censoring_rate~Beta) + 
  theme(strip.text.x = element_text(size = 11, color = "black"),
        strip.text.y = element_text(size = 11, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 12, color = "black"))

meta_specificity_plot <- Results_all_feature_selection_metabolome_summary %>% 
  filter(Censoring_rate %in% c("0.1", "0.3", "0.5", "0.7")) %>%
  ggplot(aes(x = Noise_to_signal_ratio, y = Specificity, color = Method, shape = Effect_scale, group = Method_effect)) + 
  theme(legend.position = "none") +
  geom_line() + geom_point() +
  xlab("Noise") + guides(fill = "none") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  facet_grid(Censoring_rate~Beta) + 
  theme(strip.text.x = element_text(size = 11, color = "black"),
        strip.text.y = element_text(size = 11, color = "black"),
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 12, color = "black"))

measure_plot <- ggarrange(micro_sensitivity_plot, meta_sensitivity_plot, micro_specificity_plot, meta_specificity_plot, 
                       nrow = 2, ncol = 2,  labels = c("A", "B", "C", "D"), common.legend = T)

jpeg(file = "Measure_plot_ggarrange.jpeg", width = 3000, height = 2500, res = 250, quality = 1200)
measure_plot
dev.off()
