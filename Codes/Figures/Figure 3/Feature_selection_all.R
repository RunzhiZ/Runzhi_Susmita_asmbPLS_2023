Results_all_feature_selection_microbiome <- NULL
Results_all_feature_selection_metabolome <- NULL

dist <- "lognormal"
q = 200
p = 200

c = "00"
r = "00"
effect = "05"

for (c in c("00", "01", "03", "05", "06", "07")) {
  for (r in c("00", "01", "02", "05", "10")) {
    for (effect in c("05", "20")) {
      setwd(paste0("C:/OneDrive/OneDrive - University of Florida/My research/asmbPLS/PLS/Simulation Hiper asmbPLS/q.", q, ".p.", p, "/", dist, "_mean_imputation/Results_V4"))
      eval(parse(text = paste0("load(\"Results_c", c, "_r", r, "_e", effect, "_", dist, ".RData\")")))
      for (n in c(1:3, 5, 6)) {
        #### feature selection
        temp_1 <- as.data.frame(matrix(NA, nrow = 400, ncol = 11))
        temp_2 <- as.data.frame(matrix(NA, nrow = 400, ncol = 11))
        colnames(temp_1) <- colnames(temp_2) <- c("Distribution", "Censoring_rate", "Noise_to_signal_ratio", "Effect_scale",
                            "Beta", "Method", "True", "TP", "TN", "FP", "FN")
        c_input <- c("0", "0.1", "0.3", "0.5", "0.6", "0.7")[which(c("00", "01", "03", "05", "06", "07") %in% c)]
        r_input <- c("0", "0.1", "0.2", "0.5", "1")[which(c("00", "01", "02", "05", "10") %in% r)]
        effect_input <- c("0.5", "2")[which(c("05", "20") %in% effect)]
        
        temp_1$Censoring_rate <- c_input
        temp_2$Censoring_rate <- c_input
        temp_1$Noise_to_signal_ratio <- r_input
        temp_2$Noise_to_signal_ratio <- r_input
        temp_1$Effect_scale <- effect_input
        temp_2$Effect_scale <- effect_input
        temp_1$Beta <- paste0("beta ", n)
        temp_2$Beta <- paste0("beta ", n)
        temp_1$Method <- rep(c("asmbPLS", "IPF", "SGL", "Priority"), each = 100)
        temp_2$Method <- rep(c("asmbPLS", "IPF", "SGL", "Priority"), each = 100)
        
        if(any(c(1, 5, 6) %in% n)){
          True_number_1 <- c(5, 10, 5)[which(c(1, 5, 6) %in% n)]
          True_number_2 <- c(5, 5, 10)[which(c(1, 5, 6) %in% n)]
          temp_1$True <- True_number_1
          temp_2$True <- True_number_2
          temp_asmbPLS_1_TP <- temp_asmbPLS_1_TN <- temp_asmbPLS_2_TP <- temp_asmbPLS_2_TN <- NULL
          temp_IPF_1_TP <- temp_IPF_1_TN <- temp_IPF_2_TP <- temp_IPF_2_TN <- NULL
          temp_SGL_1_TP <- temp_SGL_1_TN <- temp_SGL_2_TP <- temp_SGL_2_TN <- NULL
          temp_Priority_1_TP <- temp_Priority_1_TN <- temp_Priority_2_TP <- temp_Priority_2_TN <- NULL
          for(i in 1:100) {
            ## True positive
            temp_asmbPLS_1_TP <- c(temp_asmbPLS_1_TP, sum(results_feature_selection_list[[n]][[i]][[1]][1:True_number_1, "asmbPLS_0.7"]!=0))
            temp_asmbPLS_2_TP <- c(temp_asmbPLS_2_TP, sum(results_feature_selection_list[[n]][[i]][[2]][1:True_number_2, "asmbPLS_0.7"]!=0))
            temp_IPF_1_TP <- c(temp_IPF_1_TP, sum(results_feature_selection_list[[n]][[i]][[1]][1:True_number_1, "IPF"]!=0))
            temp_IPF_2_TP <- c(temp_IPF_2_TP, sum(results_feature_selection_list[[n]][[i]][[2]][1:True_number_2, "IPF"]!=0))
            temp_SGL_1_TP <- c(temp_SGL_1_TP, sum(results_feature_selection_list[[n]][[i]][[1]][1:True_number_1, "SGL"]!=0))
            temp_SGL_2_TP <- c(temp_SGL_2_TP, sum(results_feature_selection_list[[n]][[i]][[2]][1:True_number_2, "SGL"]!=0))
            temp_Priority_1_TP <- c(temp_Priority_1_TP, sum(results_feature_selection_list[[n]][[i]][[1]][1:True_number_1, "priority_lasso"]!=0))
            temp_Priority_2_TP <- c(temp_Priority_2_TP, sum(results_feature_selection_list[[n]][[i]][[2]][1:True_number_2, "priority_lasso"]!=0))
            ## True Negative
            temp_asmbPLS_1_TN <- c(temp_asmbPLS_1_TN, sum(which(results_feature_selection_list[[n]][[i]][[1]][, "asmbPLS_0.7"]==0) %in% (True_number_1+1):q))
            temp_asmbPLS_2_TN <- c(temp_asmbPLS_2_TN, sum(which(results_feature_selection_list[[n]][[i]][[1]][, "asmbPLS_0.7"]==0) %in% (True_number_2+1):p))
            temp_IPF_1_TN <- c(temp_IPF_1_TN, sum(which(results_feature_selection_list[[n]][[i]][[1]][, "IPF"]==0) %in% (True_number_1+1):q))
            temp_IPF_2_TN <- c(temp_IPF_2_TN, sum(which(results_feature_selection_list[[n]][[i]][[1]][, "IPF"]==0) %in% (True_number_2+1):p))
            temp_SGL_1_TN <- c(temp_SGL_1_TN, sum(which(results_feature_selection_list[[n]][[i]][[1]][, "SGL"]==0) %in% (True_number_1+1):q))
            temp_SGL_2_TN <- c(temp_SGL_2_TN, sum(which(results_feature_selection_list[[n]][[i]][[1]][, "SGL"]==0) %in% (True_number_2+1):p))
            temp_Priority_1_TN <- c(temp_Priority_1_TN, sum(which(results_feature_selection_list[[n]][[i]][[1]][, "priority_lasso"]==0) %in% (True_number_1+1):q))
            temp_Priority_2_TN <- c(temp_Priority_2_TN, sum(which(results_feature_selection_list[[n]][[i]][[1]][, "priority_lasso"]==0) %in% (True_number_2+1):p))
          }
          ## True positive
          temp_1$TP <- c(temp_asmbPLS_1_TP, temp_IPF_1_TP, temp_SGL_1_TP, temp_Priority_1_TP)
          temp_2$TP <- c(temp_asmbPLS_2_TP, temp_IPF_2_TP, temp_SGL_2_TP, temp_Priority_2_TP)
          ## True Negative
          temp_1$TN <- c(temp_asmbPLS_1_TN, temp_IPF_1_TN, temp_SGL_1_TN, temp_Priority_1_TN)
          temp_2$TN <- c(temp_asmbPLS_2_TN, temp_IPF_2_TN, temp_SGL_2_TN, temp_Priority_2_TN)
          ## False positive
          temp_1$FP <- (q - temp_1$True) - temp_1$TN
          temp_2$FP <- (p - temp_2$True) - temp_2$TN
          ## False negative
          temp_1$FN <- temp_1$True - temp_1$TP
          temp_2$FN <- temp_2$True - temp_2$TP
        }
        if(any(c(2,3) %in% n)){
          temp_true_1 <- temp_true_2 <- NULL
          temp_asmbPLS_1_TP <- temp_asmbPLS_1_TN <- temp_asmbPLS_2_TP <- temp_asmbPLS_2_TN <- NULL
          temp_IPF_1_TP <- temp_IPF_1_TN <- temp_IPF_2_TP <- temp_IPF_2_TN <- NULL
          temp_SGL_1_TP <- temp_SGL_1_TN <- temp_SGL_2_TP <- temp_SGL_2_TN <- NULL
          temp_Priority_1_TP <- temp_Priority_1_TN <- temp_Priority_2_TP <- temp_Priority_2_TN <- NULL
          for(i in 1:100) {
            temp_true_1 <- c(temp_true_1, length(which(results_feature_selection_list[[n]][[i]][[1]][, "p_value_adjusted"] <= 0.05)))
            temp_true_2 <- c(temp_true_2, length(which(results_feature_selection_list[[n]][[i]][[2]][, "p_value_adjusted"] <= 0.05)))
            index_true_1 <- which(results_feature_selection_list[[n]][[i]][[1]][, "p_value_adjusted"] <= 0.05)
            index_true_2 <- which(results_feature_selection_list[[n]][[i]][[2]][, "p_value_adjusted"] <= 0.05)
            ## True positive
            temp_asmbPLS_1_TP <- c(temp_asmbPLS_1_TP, sum(results_feature_selection_list[[n]][[i]][[1]][index_true_1, "asmbPLS_0.7"]!=0))
            temp_asmbPLS_2_TP <- c(temp_asmbPLS_2_TP, sum(results_feature_selection_list[[n]][[i]][[2]][index_true_2, "asmbPLS_0.7"]!=0))
            temp_IPF_1_TP <-  c(temp_IPF_1_TP, sum(results_feature_selection_list[[n]][[i]][[1]][index_true_1, "IPF"]!=0))
            temp_IPF_2_TP <-  c(temp_IPF_2_TP, sum(results_feature_selection_list[[n]][[i]][[2]][index_true_2, "IPF"]!=0))
            temp_SGL_1_TP <-  c(temp_SGL_1_TP, sum(results_feature_selection_list[[n]][[i]][[1]][index_true_1, "SGL"]!=0))
            temp_SGL_2_TP <-  c(temp_SGL_2_TP, sum(results_feature_selection_list[[n]][[i]][[2]][index_true_2, "SGL"]!=0))
            temp_Priority_1_TP <-  c(temp_Priority_1_TP, sum(results_feature_selection_list[[n]][[i]][[1]][index_true_1, "priority_lasso"]!=0))
            temp_Priority_2_TP <-  c(temp_Priority_2_TP, sum(results_feature_selection_list[[n]][[i]][[2]][index_true_2, "priority_lasso"]!=0))
            ## True negative
            temp_asmbPLS_1_TN <- c(temp_asmbPLS_1_TN, sum(!(which(results_feature_selection_list[[n]][[i]][[1]][, "asmbPLS_0.7"]==0) %in% index_true_1)))
            temp_asmbPLS_2_TN <- c(temp_asmbPLS_2_TN, sum(!(which(results_feature_selection_list[[n]][[i]][[2]][, "asmbPLS_0.7"]==0) %in% index_true_2)))
            temp_IPF_1_TN <-  c(temp_IPF_1_TN, sum(!(which(results_feature_selection_list[[n]][[i]][[1]][, "IPF"]==0) %in% index_true_1)))
            temp_IPF_2_TN <-  c(temp_IPF_2_TN, sum(!(which(results_feature_selection_list[[n]][[i]][[2]][, "IPF"]==0) %in% index_true_2)))
            temp_SGL_1_TN <-  c(temp_SGL_1_TN, sum(!(which(results_feature_selection_list[[n]][[i]][[1]][, "SGL"]==0) %in% index_true_1)))
            temp_SGL_2_TN <-  c(temp_SGL_2_TN, sum(!(which(results_feature_selection_list[[n]][[i]][[2]][, "SGL"]==0) %in% index_true_2)))
            temp_Priority_1_TN <-  c(temp_Priority_1_TN, sum(!(which(results_feature_selection_list[[n]][[i]][[1]][, "priority_lasso"]==0) %in% index_true_1)))
            temp_Priority_2_TN <-  c(temp_Priority_2_TN, sum(!(which(results_feature_selection_list[[n]][[i]][[2]][, "priority_lasso"]==0) %in% index_true_2)))
          }
          ## True
          temp_1$True <- temp_true_1
          temp_2$True <- temp_true_2
          ## True positive 
          temp_1$TP <- c(temp_asmbPLS_1_TP, temp_IPF_1_TP, temp_SGL_1_TP, temp_Priority_1_TP)
          temp_2$TP <- c(temp_asmbPLS_2_TP, temp_IPF_2_TP, temp_SGL_2_TP, temp_Priority_2_TP)
          ## True negative
          temp_1$TN <- c(temp_asmbPLS_1_TN, temp_IPF_1_TN, temp_SGL_1_TN, temp_Priority_1_TN)
          temp_2$TN <- c(temp_asmbPLS_2_TN, temp_IPF_2_TN, temp_SGL_2_TN, temp_Priority_2_TN)
          ## False positive
          temp_1$FP <- (q - temp_1$True) - temp_1$TN
          temp_2$FP <- (p - temp_2$True) - temp_2$TN
          ## False negative
          temp_1$FN <- temp_1$True - temp_1$TP
          temp_2$FN <- temp_2$True - temp_2$TP
        }
        Results_all_feature_selection_microbiome <- rbind(Results_all_feature_selection_microbiome, temp_1)
        Results_all_feature_selection_metabolome <- rbind(Results_all_feature_selection_metabolome, temp_2)
      }
    }
  }
}

path <- paste0("C:/OneDrive/OneDrive - University of Florida/My research/asmbPLS/PLS/Simulation Hiper asmbPLS/q.", q, ".p.", p, "/", dist, "_mean_imputation/AUC_plots")
file_name <- paste0("Results_feature_selection_", "q.", q, ".p.", p, ".RData")
setwd(path)
save(list = c(
  "Results_all_feature_selection_microbiome",
  "Results_all_feature_selection_metabolome"
),
file = file_name)
