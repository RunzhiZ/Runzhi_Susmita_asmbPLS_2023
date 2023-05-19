Results_all_for_plot <- NULL

dist <- "lognormal"
q = 200
p = 200

good_stop <- function(data, p) {
  i = 1
  while((data[i]*(1-p) > data[i+1]) & (i < 5)) {
    i = i + 1
  }
  return(i)
}


get_k <- function(data,k) {
  output <- NULL
  for (i in 1:ncol(data)) {
    output[i] <- data[k[i], i]
  }
  return(output)
}

for (c in c("00", "01", "03", "05", "06", "07")) {
  for (r in c("00", "01", "02", "05", "10")) {
    for (effect in c("05", "20")) {
      setwd(paste0("C:/OneDrive/OneDrive - University of Florida/My research/asmbPLS/PLS/Simulation Hiper asmbPLS/q.", q, ".p.", p, "/", dist, "_mean_imputation/Results_V4"))
      eval(parse(text = paste0("load(\"Results_c", c, "_r", r, "_e", effect, "_", dist, ".RData\")")))
      setwd(paste0("C:/OneDrive/OneDrive - University of Florida/My research/asmbPLS/PLS/Simulation Hiper asmbPLS/q.", q, ".p.", p, "/", dist, "_mean_imputation/Results_V4_plus"))
      eval(parse(text = paste0("load(\"Results_c", c, "_r", r, "_e", effect, "_", dist, ".RData\")")))
      for (n in 1:6) {
        temp <- as.data.frame(matrix(NA, nrow = 600, ncol = 9))
        colnames(temp) <- c("Distribution", "Censoring_rate", "Noise_to_signal_ratio", "Effect_scale",
                            "Beta", "Sigma_sqr", "Method", "MSE", "MSE_scaled")
        c_input <- c("0", "0.1", "0.3", "0.5", "0.6", "0.7")[which(c("00", "01", "03", "05", "06", "07") %in% c)]
        r_input <- c("0", "0.1", "0.2", "0.5", "1")[which(c("00", "01", "02", "05", "10") %in% r)]
        effect_input <- c("0.5", "2")[which(c("05", "20") %in% effect)]
        temp$Censoring_rate <- c_input
        temp$Noise_to_signal_ratio <- r_input
        temp$Effect_scale <- effect_input
        temp$Beta <- paste0("beta ", n)
        temp$Sigma_sqr <- rep(results_simulation_sigma_sqr_list[[n]], 6)
        temp$Method <- rep(c("asmbPLS_cv", "asmbPLS_1", "mbPLS", "BF", "IPF", "SGL"), each = 100)
        temp$MSE <- c(get_k(results_asmbPLS_0.7_0.7_predict_list[[n]], unlist(apply(results_asmbPLS_0.7_0.7_cv_list[[n]], 2, good_stop, 0.05))),
                      results_asmbPLS_0.7_0.7_predict_list[[n]][1, ], 
                      results_mbPLS_predict_list[[n]][1, ],
                      results_BF_predict_list[[n]],
                      results_IPF_predict_list[[n]],
                      results_SGL_predict_list[[n]])
        temp$MSE_scaled <- temp$MSE/temp$Sigma_sqr
        Results_all_for_plot <- rbind(Results_all_for_plot, temp)
      }
    }
  }
}
Results_all_for_plot$Method <- factor(Results_all_for_plot$Method, levels=c("asmbPLS_cv", "asmbPLS_1", "mbPLS", "BF", "IPF", "SGL"))

path <- paste0("C:/OneDrive/OneDrive - University of Florida/My research/asmbPLS/PLS/Simulation Hiper asmbPLS/q.", q, ".p.", p, "/", dist, "_mean_imputation/Boxplots")
file_name <- paste0("Results_for_plot_", "q.", q, ".p.", p, ".RData")
setwd(path)
save(list = c(
  "Results_all_for_plot"
),
file = file_name)
