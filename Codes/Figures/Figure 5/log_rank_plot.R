###############################################################################
# load R package and function
###############################################################################
library(prioritylasso)
library(glmnet)
library(blockForest)
library(ipflasso)
library(compositions)
library(Rcpp)
library(RcppArmadillo)
library(SGL)
library(readxl)
library(dplyr)
library(survminer)
library(survival)
library(asmbPLS)

###############################################################################
# read real data
###############################################################################
setwd("C:/OneDrive/OneDrive - University of Florida/My research/asmbPLS/Dietary fiber and probiotics influence the gut microbiome and melanoma immunotherapy response")
metadata <- read_excel("PD1_Wargo_Human_WGS_Relabund_and_metadata_light_filtering.xlsx", sheet = "Metadata")
microbiome_count_t <- read_excel("PD1_Wargo_Human_WGS_Relabund_and_metadata_light_filtering.xlsx", sheet = "LKT_PPM")
microbiome_taxa <- read_excel("PD1_Wargo_Human_WGS_Relabund_and_metadata_light_filtering.xlsx", sheet = "LKT_featuretable")
protein_count_t <- read_excel("PD1_Wargo_Human_WGS_Relabund_and_metadata_light_filtering.xlsx", sheet = "Product_PPM")
protein_taxa <- read_excel("PD1_Wargo_Human_WGS_Relabund_and_metadata_light_filtering.xlsx", sheet = "Product_featuretable")

metadata[(metadata == "N_A")] <- NA
## extract samples with survival data
survival_index <- complete.cases(metadata[,c("pfsevent", "pfs_d")])
metadata_survival <- metadata[survival_index, ]

###############################################################################
# pre-processing:
# 1) pre-processing for omics data
# 2) pre-processing for meta data
###############################################################################
#### pre-processing for omics data ####
## microbiome data
microbiome_count <- t(microbiome_count_t)
colnames(microbiome_count) <- microbiome_count[1,]
microbiome_count <- microbiome_count[-1,]
microbiome_count <- apply(microbiome_count, 2, as.numeric)
# pseudo value imputation 0 -> 0.5
microbiome_count[microbiome_count == 0] <- 0.5
# count -> relative abundance
microbiome_relative_abundance <- t(apply(microbiome_count, 1, function(data){data/sum(data)}))
microbiome_relative_abundance_survival <- microbiome_relative_abundance[survival_index, ]
dim(microbiome_relative_abundance_survival)
## protein data
protein_count <- t(protein_count_t)
colnames(protein_count) <- protein_count[1,]
protein_count <- protein_count[-1,]
protein_count <- apply(protein_count, 2, as.numeric)
# pseudo value imputation 0 -> 0.5
protein_count[protein_count == 0] <- 0.5
protein_count_survival <- protein_count[survival_index, ]
dim(protein_count_survival)

#### pre-processing for meta data ####
covariates <- metadata_survival[, c("Age", "Sex", "metformin", "steroids", "statins", "PPI", "beta_blocker", "other_anti_hypertensive",
                                    "BMI", "Treatment_naive", "LDH", "adv_substage", "response", "treatment", "probiotics", "antibiotics", "fiber_cat")]
## Age
covariates$Age[79] <- 90
covariates$Age <- as.numeric(covariates$Age)
## Sex
covariates$Sex[covariates$Sex == "F"] <- 0
covariates$Sex[covariates$Sex == "M"] <- 1
covariates$Sex <- as.numeric(covariates$Sex)
## BMI
covariates$BMI <- as.numeric(covariates$BMI)
## metformin
covariates$metformin[covariates$metformin == "No"] <- 0
covariates$metformin[covariates$metformin == "Yes"] <- 1
covariates$metformin <- as.numeric(covariates$metformin)
## steroids
covariates$steroids[covariates$steroids == "No"] <- 0
covariates$steroids[covariates$steroids == "Yes"] <- 1
covariates$steroids <- as.numeric(covariates$steroids)
## statins
covariates$statins[covariates$statins == "No"] <- 0
covariates$statins[covariates$statins == "Yes"] <- 1
covariates$statins <- as.numeric(covariates$statins)
## PPI
covariates$PPI[covariates$PPI == "No"] <- 0
covariates$PPI[covariates$PPI == "Yes"] <- 1
covariates$PPI <- as.numeric(covariates$PPI)
## beta_blocker
covariates$beta_blocker[covariates$beta_blocker == "No"] <- 0
covariates$beta_blocker[covariates$beta_blocker == "Yes"] <- 1
covariates$beta_blocker <- as.numeric(covariates$beta_blocker)
## other_anti_hypertensive
covariates$other_anti_hypertensive[covariates$other_anti_hypertensive == "No"] <- 0
covariates$other_anti_hypertensive[covariates$other_anti_hypertensive == "Yes"] <- 1
covariates$other_anti_hypertensive <- as.numeric(covariates$other_anti_hypertensive)
## Treatment_naive
covariates$Treatment_naive[covariates$Treatment_naive == "No"] <- 0
covariates$Treatment_naive[covariates$Treatment_naive == "Yes"] <- 1
covariates$Treatment_naive <- as.numeric(covariates$Treatment_naive)
## LDH
covariates$LDH[covariates$LDH == "No"] <- 0
covariates$LDH[covariates$LDH == "Yes"] <- 1
covariates$LDH <- as.numeric(covariates$LDH)
## adv_substage
covariates$adv_substage[covariates$adv_substage == "Stage_M1C"] <- 0
covariates$adv_substage[covariates$adv_substage == "Stage_M1D"] <- 1
covariates$adv_substage <- as.numeric(covariates$adv_substage)
## response
covariates$response[covariates$response == "Non_Responder"] <- 0
covariates$response[covariates$response == "Responder"] <- 1
covariates$response <- as.numeric(covariates$response)
## treatment
covariates$treatment[covariates$treatment == "all ICB"] <- 0
covariates$treatment[covariates$treatment == "anti-PD1"] <- 1
covariates$treatment <- as.numeric(covariates$treatment)
## probiotics
covariates$probiotics[covariates$probiotics == "No"] <- 0
covariates$probiotics[covariates$probiotics == "Yes"] <- 1
covariates$probiotics <- as.numeric(covariates$probiotics)
## antibiotics
covariates$antibiotics[covariates$antibiotics == "No"] <- 0
covariates$antibiotics[covariates$antibiotics == "Yes"] <- 1
covariates$antibiotics <- as.numeric(covariates$antibiotics)
## fiber_cat
covariates$fiber_cat[covariates$fiber_cat == "Insufficient_intake"] <- 0
covariates$fiber_cat[covariates$fiber_cat == "Sufficient_intake"] <- 1
covariates$fiber_cat <- as.numeric(covariates$fiber_cat)

covatiates_index <- complete.cases(covariates)

###############################################################################
# quantile combination setting
###############################################################################
quantile_1 <- c(0.9, 0.925, 0.95, 0.975, 0.99, 0.999)
quantile_2 <- c(0.997, 0.9985, 0.9993, 0.9999)
quantile_3 <- c(0, 0.3, 0.5, 0.7, 0.8, 0.9, 0.99)
quantilelist <- list(quantile_1, quantile_2)
quantilelist_covariates <- list(quantile_1, quantile_2, quantile_3)
quantile_table <- quantileComb(quantilelist)
quantile_table_covariates <- quantileComb(quantilelist_covariates)

###############################################################################
# prepare X and Y matrix for model fitting
###############################################################################
#### without covariates ####
X_dim <- c(ncol(microbiome_relative_abundance_survival), ncol(protein_count_survival))
X <- cbind(microbiome_relative_abundance_survival, protein_count_survival)
## imputation method
survival_time <- metadata_survival[, c("pfs_d", "pfsevent")]
survival_time <- survival_time %>% mutate_at(c("pfs_d", "pfsevent"), as.numeric)
survival_time$pfs_d_log <- log(survival_time$pfs_d)
survival_time_imp <- meanimp(survival_time[,c("pfs_d_log", "pfsevent")])$imputed_table
Y <- as.matrix(survival_time_imp[, 4])
Y_indicator <- survival_time_imp[, 2]

#### with covariates ####
X_dim <- c(ncol(microbiome_relative_abundance_survival), ncol(protein_count_survival), ncol(covariates))
X <- as.matrix(cbind(microbiome_relative_abundance_survival, protein_count_survival, covariates))
X <- X[covatiates_index, ]
## imputation method
survival_time <- metadata_survival[covatiates_index, c("pfs_d", "pfsevent")]
survival_time <- survival_time %>% mutate_at(c("pfs_d", "pfsevent"), as.numeric)
survival_time$pfs_d_log <- log(survival_time$pfs_d)
survival_time_imp <- meanimp(survival_time[,c("pfs_d_log", "pfsevent")])$imputed_table
Y <- as.matrix(survival_time_imp[, 4])
Y_indicator <- survival_time_imp[, 2]

## clr transformation for X ##
X_clr <- X
X_clr[, 1:X_dim[1]] <- clr(X[, 1:X_dim[1]])
X_clr[, (X_dim[1]+1):(X_dim[1]+X_dim[2])] <- clr(X[, (X_dim[1]+1):(X_dim[1]+X_dim[2])])

#### compare function: compare predicted and true Y ####
result_compare <- function (prediction, true) {
  result_abs = sum(abs(prediction-true))/length(true)
  result_sq = sum((prediction-true)^2)/length(true)
  return(list(result_sq = result_sq,
              result_abs = result_abs))
}

observed_index <- which(Y_indicator == 1)

#### asmbPLS ####
PLS_term_selected <- 3
k = 5
## MSE_fit ##
asmbPLS_results_cv <- asmbPLS.cv(X.matrix = X_clr, 
                                 Y.matrix = Y, 
                                 PLS.comp = PLS_term_selected,
                                 X.dim = X_dim, 
                                 quantile.comb.table = quantile_table_covariates[,1:length(X_dim)], 
                                 Y.indicator = Y_indicator ,
                                 k = k)
PLS_table <- asmbPLS_results_cv$quantile_table_CV[,1:length(X_dim)]
asmbPLS_results <- asmbPLS.fit(X.matrix = X_clr,
                               Y.matrix = Y, 
                               PLS.comp = 3,
                               X.dim = X_dim, 
                               quantile.comb = PLS_table[1:PLS_term_selected,])
Y_fit_asmbPLS <- asmbPLS.predict(fit.results = asmbPLS_results,
                                 X.matrix.new = X_clr, 
                                 PLS.comp = 1)
result_compare(Y_fit_asmbPLS$Y_pred[observed_index], Y[observed_index])[[1]]

## super score for survival ##
survival_time_imp_ss <- cbind(survival_time_imp, asmbPLS_results$X_super_score[,1])
colnames(survival_time_imp_ss)[ncol(survival_time_imp_ss)] <- "super_score"
cutoff_super_score <- surv_cutpoint(
  survival_time_imp_ss,
  time = "Survival_time",
  event = "Event_indicator",
  variables = c("super_score")
)
survival_time_imp_ss$group <- 0
survival_time_imp_ss$group[survival_time_imp_ss$super_score > cutoff_super_score$cutpoint[, 1]] <- 1
survival_time$group <- survival_time_imp_ss$group

## log-rank test ##
setwd("C:/OneDrive/OneDrive - University of Florida/My research/asmbPLS/Dietary fiber and probiotics influence the gut microbiome and melanoma immunotherapy response")
logr_fit <- survfit(Surv(pfs_d, pfsevent) ~ group, data = survival_time)
logrank_plot <- ggsurvplot(logr_fit,
                           pval = TRUE, conf.int = TRUE,
                           risk.table = F, # Add risk table
                           risk.table.col = "strata", # Change risk table color by groups
                           linetype = "strata", # Change line type by groups
                           surv.median.line = "hv", # Specify median survival
                           ggtheme = theme_bw(), # Change ggplot2 theme
                           palette = c("#E7B800", "#2E9FDF"))
jpeg(file = paste0("logrank_plot.jpeg"), width = 2000, height = 1500, res = 250, quality = 1200)
logrank_plot
dev.off()
