library(ggplot2)
library(dplyr)
library(pROC)

dist <- "lognormal"
q = 200
p = 200
path <- paste0("C:/OneDrive/OneDrive - University of Florida/My research/asmbPLS/PLS/Simulation Hiper asmbPLS/q.", q, ".p.", p, "/", dist, "_mean_imputation/AUC_plots")
setwd(path)
file_name <- paste0("Results_feature_selection_", "q.", q, ".p.", p, ".RData")
load(file_name)

Results_all_feature_selection_microbiome_summary <- Results_all_feature_selection_microbiome %>% distinct(Censoring_rate, Noise_to_signal_ratio, Effect_scale, Beta, Method, .keep_all = T)
Results_all_feature_selection_metabolome_summary <- Results_all_feature_selection_metabolome %>% distinct(Censoring_rate, Noise_to_signal_ratio, Effect_scale, Beta, Method, .keep_all = T)

for(i in 1:nrow(Results_all_feature_selection_microbiome_summary)) {
  index_1 <- which((Results_all_feature_selection_microbiome$Censoring_rate == Results_all_feature_selection_microbiome_summary$Censoring_rate[i])
                   &(Results_all_feature_selection_microbiome$Noise_to_signal_ratio == Results_all_feature_selection_microbiome_summary$Noise_to_signal_ratio[i])
                   &(Results_all_feature_selection_microbiome$Effect_scale == Results_all_feature_selection_microbiome_summary$Effect_scale[i])
                   &(Results_all_feature_selection_microbiome$Beta == Results_all_feature_selection_microbiome_summary$Beta[i])
                   &(Results_all_feature_selection_microbiome$Method == Results_all_feature_selection_microbiome_summary$Method[i])
  )
  index_2 <- which((Results_all_feature_selection_metabolome$Censoring_rate == Results_all_feature_selection_metabolome_summary$Censoring_rate[i])
                   &(Results_all_feature_selection_metabolome$Noise_to_signal_ratio == Results_all_feature_selection_metabolome_summary$Noise_to_signal_ratio[i])
                   &(Results_all_feature_selection_metabolome$Effect_scale == Results_all_feature_selection_metabolome_summary$Effect_scale[i])
                   &(Results_all_feature_selection_metabolome$Beta == Results_all_feature_selection_metabolome_summary$Beta[i])
                   &(Results_all_feature_selection_metabolome$Method == Results_all_feature_selection_metabolome_summary$Method[i])
  )
  temp_1 <- Results_all_feature_selection_microbiome[index_1, ]
  temp_2 <- Results_all_feature_selection_metabolome[index_2, ]
  Results_all_feature_selection_microbiome_summary[i, c("True", "TP", "TN", "FP", "FN")] <- apply(Results_all_feature_selection_microbiome[index_1, c("True", "TP", "TN", "FP", "FN")], 2, sum)
  Results_all_feature_selection_metabolome_summary[i, c("True", "TP", "TN", "FP", "FN")] <- apply(Results_all_feature_selection_metabolome[index_1, c("True", "TP", "TN", "FP", "FN")], 2, sum)
}

Results_all_feature_selection_microbiome_summary$Sensitivity <- Results_all_feature_selection_microbiome_summary$TP/(Results_all_feature_selection_microbiome_summary$TP + Results_all_feature_selection_microbiome_summary$FN)
Results_all_feature_selection_microbiome_summary$Specificity <- Results_all_feature_selection_microbiome_summary$TN/(Results_all_feature_selection_microbiome_summary$TN + Results_all_feature_selection_microbiome_summary$FP)

Results_all_feature_selection_metabolome_summary$Sensitivity <- Results_all_feature_selection_metabolome_summary$TP/(Results_all_feature_selection_metabolome_summary$TP + Results_all_feature_selection_metabolome_summary$FN)
Results_all_feature_selection_metabolome_summary$Specificity <- Results_all_feature_selection_metabolome_summary$TN/(Results_all_feature_selection_metabolome_summary$TN + Results_all_feature_selection_metabolome_summary$FP)

## AUC
Results_all_feature_selection_microbiome_summary$AUC <- NA
Results_all_feature_selection_metabolome_summary$AUC <- NA
for(i in 1:nrow(Results_all_feature_selection_microbiome_summary)) {
  row_1 <- Results_all_feature_selection_microbiome_summary[i, ]
  row_2 <- Results_all_feature_selection_metabolome_summary[i, ]
  temp_1 <- cbind(c(rep(1, row_1$TP), rep(0, row_1$TN), rep(0, row_1$FP), rep(1, row_1$FN)), 
                  c(rep(1, row_1$TP), rep(0, row_1$TN), rep(1, row_1$FP), rep(0, row_1$FN)))
  temp_2 <- cbind(c(rep(1, row_2$TP), rep(0, row_2$TN), rep(0, row_2$FP), rep(1, row_2$FN)), 
                  c(rep(1, row_2$TP), rep(0, row_2$TN), rep(1, row_2$FP), rep(0, row_2$FN)))
  Results_all_feature_selection_microbiome_summary$AUC[i] <- roc(temp_1[, 1], temp_1[, 2], quiet = TRUE)$auc
  Results_all_feature_selection_metabolome_summary$AUC[i] <- roc(temp_2[, 1], temp_2[, 2], quiet = TRUE)$auc
}

# p_sensitivity <- Results_all_feature_selection_microbiome_summary %>% 
#   filter(Censoring_rate == 0.3 & Effect_scale == 0.5) %>%
#   mutate(Type = ifelse(Method == "asmbPLS", "asmbPLS", "Others")) %>%
#   ggplot(aes(x = Method, y = Sensitivity, fill = Type)) + 
#   #geom_jitter(width=0.1) +
#   scale_fill_manual(values=c("grey", "white")) +
#   theme(legend.position = "none") +
#   geom_boxplot() + xlab("") + 
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust=1)) +
#   facet_grid(Noise_to_signal_ratio~Beta)
# p_specificity <- Results_all_feature_selection_microbiome_summary %>% 
#   filter(Censoring_rate == 0.5 & Effect_scale == 2) %>%
#   mutate(Type = ifelse(Method == "asmbPLS", "asmbPLS", "Others")) %>%
#   ggplot(aes(x = Method, y = Specificity, fill = Type)) + 
#   #geom_jitter(width=0.1) +
#   scale_fill_manual(values=c("grey", "white")) +
#   theme(legend.position = "none") +
#   geom_boxplot() + xlab("") + 
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust=1)) +
#   facet_grid(Noise_to_signal_ratio~Beta)

save.image("Feature_selection_summary.RData")

AUC_plot_generate <- function(data, c, e) {
  AUC_plot <- data %>% 
    filter(Censoring_rate == c & Effect_scale == e) %>%
    mutate(Methods = ifelse(Method == "asmbPLS", "asmbPLS", "Others")) %>%
    ggplot(aes(x = Method, y = AUC, fill = Methods)) + 
    #geom_jitter(width=0.1) +
    scale_fill_manual(values=c("black", "grey")) +
    theme(legend.position = "none") +
    geom_bar(stat = 'identity') + xlab("") + guides(fill = "none") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust=1)) +
    facet_grid(Noise_to_signal_ratio~Beta)
  return(AUC_plot)
}

#### microbiome, 0.5 ####
path <- paste0("C:/OneDrive/OneDrive - University of Florida/My research/asmbPLS/PLS/Simulation Hiper asmbPLS/q.", q, ".p.", p, "/", dist, "_mean_imputation/AUC_plots/microbiome/0.5")
setwd(path)
c = 0.1
e = 0.5
AUC_plot <- AUC_plot_generate(Results_all_feature_selection_microbiome_summary, c, e)
jpeg(file = paste0("p_AUC_", c, "_", e, "_micro.jpeg"), width = 1200, height = 800, res = 120, quality = 200)
AUC_plot
dev.off()

c = 0.3
e = 0.5
AUC_plot <- AUC_plot_generate(Results_all_feature_selection_microbiome_summary, c, e)
jpeg(file = paste0("p_AUC_", c, "_", e, "_micro.jpeg"), width = 1200, height = 800, res = 120, quality = 200)
AUC_plot
dev.off()

c = 0.5
e = 0.5
AUC_plot <- AUC_plot_generate(Results_all_feature_selection_microbiome_summary, c, e)
jpeg(file = paste0("p_AUC_", c, "_", e, "_micro.jpeg"), width = 1200, height = 800, res = 120, quality = 200)
AUC_plot
dev.off()

c = 0.7
e = 0.5
AUC_plot <- AUC_plot_generate(Results_all_feature_selection_microbiome_summary, c, e)
jpeg(file = paste0("p_AUC_", c, "_", e, "_micro.jpeg"), width = 1200, height = 800, res = 120, quality = 200)
AUC_plot
dev.off()

#### microbiome, 2 ####
path <- paste0("C:/OneDrive/OneDrive - University of Florida/My research/asmbPLS/PLS/Simulation Hiper asmbPLS/q.", q, ".p.", p, "/", dist, "_mean_imputation/AUC_plots/microbiome/2")
setwd(path)
c = 0.1
e = 2
AUC_plot <- AUC_plot_generate(Results_all_feature_selection_microbiome_summary, c, e)
jpeg(file = paste0("p_AUC_", c, "_", e, "_micro.jpeg"), width = 1200, height = 800, res = 120, quality = 200)
AUC_plot
dev.off()

c = 0.3
e = 2
AUC_plot <- AUC_plot_generate(Results_all_feature_selection_microbiome_summary, c, e)
jpeg(file = paste0("p_AUC_", c, "_", e, "_micro.jpeg"), width = 1200, height = 800, res = 120, quality = 200)
AUC_plot
dev.off()

c = 0.5
e = 2
AUC_plot <- AUC_plot_generate(Results_all_feature_selection_microbiome_summary, c, e)
jpeg(file = paste0("p_AUC_", c, "_", e, "_micro.jpeg"), width = 1200, height = 800, res = 120, quality = 200)
AUC_plot
dev.off()

c = 0.7
e = 2
AUC_plot <- AUC_plot_generate(Results_all_feature_selection_microbiome_summary, c, e)
jpeg(file = paste0("p_AUC_", c, "_", e, "_micro.jpeg"), width = 1200, height = 800, res = 120, quality = 200)
AUC_plot
dev.off()

#### metabolome, 0.5 ####
path <- paste0("C:/OneDrive/OneDrive - University of Florida/My research/asmbPLS/PLS/Simulation Hiper asmbPLS/q.", q, ".p.", p, "/", dist, "_mean_imputation/AUC_plots/metabolome/0.5")
setwd(path)
c = 0.1
e = 0.5
AUC_plot <- AUC_plot_generate(Results_all_feature_selection_metabolome_summary, c, e)
jpeg(file = paste0("p_AUC_", c, "_", e, "_meta.jpeg"), width = 1200, height = 800, res = 120, quality = 200)
AUC_plot
dev.off()

c = 0.3
e = 0.5
AUC_plot <- AUC_plot_generate(Results_all_feature_selection_metabolome_summary, c, e)
jpeg(file = paste0("p_AUC_", c, "_", e, "_meta.jpeg"), width = 1200, height = 800, res = 120, quality = 200)
AUC_plot
dev.off()

c = 0.5
e = 0.5
AUC_plot <- AUC_plot_generate(Results_all_feature_selection_metabolome_summary, c, e)
jpeg(file = paste0("p_AUC_", c, "_", e, "_meta.jpeg"), width = 1200, height = 800, res = 120, quality = 200)
AUC_plot
dev.off()

c = 0.7
e = 0.5
AUC_plot <- AUC_plot_generate(Results_all_feature_selection_metabolome_summary, c, e)
jpeg(file = paste0("p_AUC_", c, "_", e, "_meta.jpeg"), width = 1200, height = 800, res = 120, quality = 200)
AUC_plot
dev.off()

#### metabolome, 2 ####
path <- paste0("C:/OneDrive/OneDrive - University of Florida/My research/asmbPLS/PLS/Simulation Hiper asmbPLS/q.", q, ".p.", p, "/", dist, "_mean_imputation/AUC_plots/metabolome/2")
setwd(path)
c = 0.1
e = 2
AUC_plot <- AUC_plot_generate(Results_all_feature_selection_metabolome_summary, c, e)
jpeg(file = paste0("p_AUC_", c, "_", e, "_meta.jpeg"), width = 1200, height = 800, res = 120, quality = 200)
AUC_plot
dev.off()

c = 0.3
e = 2
AUC_plot <- AUC_plot_generate(Results_all_feature_selection_metabolome_summary, c, e)
jpeg(file = paste0("p_AUC_", c, "_", e, "_meta.jpeg"), width = 1200, height = 800, res = 120, quality = 200)
AUC_plot
dev.off()

c = 0.5
e = 2
AUC_plot <- AUC_plot_generate(Results_all_feature_selection_metabolome_summary, c, e)
jpeg(file = paste0("p_AUC_", c, "_", e, "_meta.jpeg"), width = 1200, height = 800, res = 120, quality = 200)
AUC_plot
dev.off()

c = 0.7
e = 2
AUC_plot <- AUC_plot_generate(Results_all_feature_selection_metabolome_summary, c, e)
jpeg(file = paste0("p_AUC_", c, "_", e, "_meta.jpeg"), width = 1200, height = 800, res = 120, quality = 200)
AUC_plot
dev.off()