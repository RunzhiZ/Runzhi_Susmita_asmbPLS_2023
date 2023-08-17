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

## blocks setting ##
blocks_vector <- rep(1:length(X_dim), times = X_dim)
blocks <- lapply(1:3, function(x) which(blocks_vector == x))
blocks_1 <- lapply(1:3, function(x) which(blocks_vector == x))
blocks_2 <- lapply(c(1, 3, 2), function(x) which(blocks_vector == x))
blocks_3 <- lapply(c(2, 1, 3), function(x) which(blocks_vector == x))
blocks_4 <- lapply(c(2, 3, 1), function(x) which(blocks_vector == x))
blocks_5 <- lapply(c(3, 1, 2), function(x) which(blocks_vector == x))
blocks_6 <- lapply(c(3, 2, 1), function(x) which(blocks_vector == x))

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

###############################################################################
# comparison between different methods:
# 1) LOO, MSE_prediction

# methods include:
# 1) asmbPLS
# 2) mbPLS
# 3) block forest
# 4) IPF-lasso
# 5) SGL
# 6) priority-lasso
###############################################################################
## CV for asmbPLS ##
observed_index <- which(Y_indicator == 1)
PLS_term_selected <- 3
k = 5
asmbPLS_results_cv <- asmbPLS.cv(X.matrix = X_clr, 
                                 Y.matrix = Y, 
                                 PLS.comp = PLS_term_selected,
                                 X.dim = X_dim, 
                                 quantile.comb.table = quantile_table_covariates[,1:length(X_dim)], 
                                 Y.indicator = Y_indicator ,
                                 k = k)
PLS_table <- asmbPLS_results_cv$quantile_table_CV[, 1:length(X_dim)]

results_LOO <- matrix(NA, nrow = length(observed_index), ncol = 11)
colnames(results_LOO) <- c("index", "true", "asmbPLS_1", "asmbPLS_2", "asmbPLS_3", 
                           "mbPLS_1", "mbPLS_2", "mbPLS_3", "BF", "IPF", "SGL")
results_LOO_MSE <- matrix(NA, nrow = length(observed_index), ncol = 10)
colnames(results_LOO_MSE) <- c("index", "asmbPLS_1", "asmbPLS_2", "asmbPLS_3", 
                               "mbPLS_1", "mbPLS_2", "mbPLS_3", "BF", "IPF", "SGL")

for(i in 1:length(observed_index)){
  print(paste0("****************************** index ", i, " ******************************"))
  #### leave-one-out ####
  i_index <- observed_index[i]
  X_test <- matrix(X_clr[i_index, ], ncol = ncol(X_clr))
  Y_test <- matrix(Y[i_index,], ncol = ncol(Y))
  colnames(X_test) <- colnames(X_clr)
  X_train <- X_clr[-i_index, ]
  Y_train <- matrix(Y[-i_index,], ncol = ncol(Y))
  Y_indicator_train <- Y_indicator[-i_index]
  
  results_LOO[i, "index"] <- i_index
  results_LOO[i, "true"] <- Y_test
  
  results_LOO_MSE[i, "index"] <- i_index
  
  #### asmbPLS ####
  PLS_term_selected <- 3
  asmbPLS_results <- asmbPLS.fit(X_train, Y_train, PLS_term_selected, X_dim, PLS_table[1:PLS_term_selected,])
  results_asmbPLS_predict_single <- NULL
  results_asmbPLS_fit_single <- NULL
  for (PLS_term in 1:PLS_term_selected) {
    Y_predict_asmbPLS <- asmbPLS.predict(asmbPLS_results, X_test, PLS_term)
    results_LOO[i, PLS_term + 2] <- Y_predict_asmbPLS$Y_pred
    results_asmbPLS_predict_single[PLS_term] <- result_compare(Y_predict_asmbPLS$Y_pred, Y_test)[[1]]
  }
  results_LOO_MSE[i, c("asmbPLS_1", "asmbPLS_2", "asmbPLS_3")] <- results_asmbPLS_predict_single[1:3]
  
  #### mbPLS ####
  mbPLS_results <- mbPLS.fit(X_train, Y_train, PLS_term_selected, X_dim)
  results_mbPLS_predict_single <- NULL
  results_mbPLS_fit_single <- NULL
  for (PLS_term in 1:PLS_term_selected) {
    Y_predict_mbPLS <- asmbPLS.predict(mbPLS_results, X_test, PLS_term)
    results_LOO[i, PLS_term + 5] <- Y_predict_mbPLS$Y_pred
    results_mbPLS_predict_single[PLS_term] <- result_compare(Y_predict_mbPLS$Y_pred, Y_test)[[1]]
  }
  
  results_LOO_MSE[i, c("mbPLS_1", "mbPLS_2", "mbPLS_3")] <- results_mbPLS_predict_single[1:3]
  
  #### Block forest ####
  tryCatch({
    Y_train_vector <- as.vector(Y_train)
    blockforobj <- blockfor(X_train, Y_train_vector, replace = TRUE, blocks = blocks, block.method = "BlockForest", seed = i)
    Y_predict_BF <- predict(blockforobj$forest, data = X_test)
    results_LOO[i, "BF"] <- Y_predict_BF
    results_LOO_MSE[i, "BF"] <- result_compare(Y_predict_BF, Y_test)[[1]]
  }, error = function(e){print(paste0("index:", i, " BF cannot predict"))})
  
  #### IPF ####
  set.seed(i)
  Y_train_vector <- as.vector(Y_train)
  Y_vector <- as.vector(Y)
  pflist <- list(c(1, 1, 1), c(2, 1, 2), c(1, 2, 2), c(3, 1, 3), c(1, 3, 3), c(4, 1, 4), c(1, 4, 4), c(5, 1, 5), c(1, 5, 5))
  ipf_fit <- cvr2.ipflasso(X = X_train, Y = Y_train_vector, family = "gaussian", type.measure = "mse",
                           blocks = blocks,
                           pflist = pflist, nfolds = 5, ncv = 10)
  Y_predict_IPF <- ipflasso.predict(ipf_fit, X_test)
  results_LOO[i, "IPF"] <- Y_predict_IPF$linpredtest
  results_LOO_MSE[i, "IPF"] <- result_compare(Y_predict_IPF$linpredtest, Y_test)[[1]]
  
  #### SGL ####
  tryCatch({
    data = list(x = X_train, y = Y_train)
    index <- c(rep(1, X_dim[1]), rep(2, X_dim[2]), rep(3, X_dim[3]))
    cvFit = cvSGL(data, index, type = "linear", nfold = 5)
    Fit = SGL(data, index, type = "linear", lambdas = cvFit$lambdas)
    Y_predict_SGL <- predictSGL(Fit, X_test, which.min(cvFit$lldiff))
    results_LOO[i, "SGL"] <- Y_predict_SGL
    results_LOO_MSE[i, "SGL"] <- result_compare(Y_predict_SGL, Y_test)[[1]]
    
  }, error = function(e){print(paste0("index:", i, "SGL cannot fit"))})
}

###############################################################################
# comparison between different methods:
# 1) fit

# methods include:
# 1) asmbPLS
# 2) mbPLS
# 3) block forest
# 4) IPF-lasso
# 5) SGL
# 6) priority-lasso
###############################################################################
fit_table <- matrix(NA, nrow = 1, ncol = 9)
colnames(fit_table) <- c("asmbPLS_1", "asmbPLS_2", "asmbPLS_3",
                         "mbPLS_1", "mbPLS_2", "mbPLS_3", "BF", "IPF", "SGL")

#### features significance table ####
results_feature_selection_list <- list()
for(i in 1:length(X_dim)) {
  results_feature_selection_list[[i]] <- matrix(NA, nrow = X_dim[i], ncol = 7)
  row.names(results_feature_selection_list[[i]]) <- colnames(X_clr)[blocks[[i]]]
  colnames(results_feature_selection_list[[i]]) <- c("p_value", "p_value_adjusted", "asmbPLS", "mbPLS", "BF", "IPF", "SGL")
}
for(i in 1:length(X_dim)) {
  for(j in 1:X_dim[i]) {
    if(sum(X_clr[,blocks[[i]][j]]!=0)) {results_feature_selection_list[[i]][j, "p_value"] <- summary(lm(Y ~ X_clr[,blocks[[i]][j]]))$coefficients[2, 4]}
  }
  results_feature_selection_list[[i]][, "p_value_adjusted"] <- p.adjust(results_feature_selection_list[[i]][, "p_value"], method = "BH")
}

#### asmbPLS ####
asmbPLS_results <- asmbPLS.fit(X_clr, Y, PLS_term_selected, X_dim, PLS_table[1:PLS_term_selected,])
for(i in 1:PLS_term_selected) {
  Y_fit_asmbPLS <- asmbPLS.predict(asmbPLS_results, X_clr, i)
  fit_table[,i] <- result_compare(Y_fit_asmbPLS$Y_pred[observed_index], Y[observed_index])[[1]]
}

## feature selection
for(i in 1:length(X_dim)) {
  results_feature_selection_list[[i]][, "asmbPLS"] <- asmbPLS_results$X_weight[[i]][,1]
}

#### mbPLS ####
mbPLS_results <- mbPLS.fit(X_clr, Y, PLS_term_selected, X_dim)
for(i in 1:PLS_term_selected) {
  Y_fit_mbPLS <- asmbPLS.predict(mbPLS_results, X_clr, i)
  fit_table[,i + 3] <- result_compare(Y_fit_mbPLS$Y_pred[observed_index], Y[observed_index])[[1]]
}

## feature selection
for(i in 1:length(X_dim)) {
  results_feature_selection_list[[i]][, "mbPLS"] <- mbPLS_results$X_weight[[i]][,1]
}

#### Block forest ####
Y_vector <- as.vector(Y)
blockforobj <- blockfor(X_clr, Y_vector, replace = TRUE, blocks = blocks, block.method = "BlockForest", importance = "impurity", seed = 1)
Y_fit_BF <- blockforobj$forest$predictions
fit_table[, "BF"] <- result_compare(Y_fit_BF[observed_index], Y[observed_index])[[1]]

## feature selection
for(i in 1:length(X_dim)) {
  results_feature_selection_list[[i]][, "BF"] <- blockforobj$forest$variable.importance[blocks[[i]]]
}

#### IPF ####
set.seed(1)
Y_vector <- as.vector(Y)
pflist <- list(c(1, 1, 1), c(2, 1, 2), c(1, 2, 2), c(3, 1, 3), c(1, 3, 3), c(4, 1, 4), c(1, 4, 4), c(5, 1, 5), c(1, 5, 5))
ipf_fit <- cvr2.ipflasso(X = X_clr, Y = Y_vector, family = "gaussian", type.measure = "mse",
                         blocks = blocks,
                         pflist = pflist, nfolds = 5, ncv = 10)
Y_fit_IPF <- ipflasso.predict(ipf_fit, X_clr)
fit_table[, "IPF"] <- result_compare(Y_fit_IPF$linpredtest[observed_index], Y[observed_index])[[1]]

## feature selection
for(i in 1:length(X_dim)) {
  results_feature_selection_list[[i]][, "IPF"] <- ipf_fit$coeff[,ipf_fit$ind.bestlambda][-1][blocks[[i]]]
}

#### SGL ####
tryCatch({
  data = list(x = X_clr, y = Y)
  index <- c(rep(1, X_dim[1]), rep(2, X_dim[2]), rep(3, X_dim[3]))
  cvFit = cvSGL(data, index, type = "linear", nfold = 5)
  Fit = SGL(data, index, type = "linear", lambdas = cvFit$lambdas)
  Y_fit_SGL <- predictSGL(Fit, X_clr, which.min(cvFit$lldiff))
  
  ## feature selection
  for(i in 1:length(X_dim)) {
    results_feature_selection_list[[i]][, "SGL"] <- Fit$beta[, which.min(cvFit$lldiff)][blocks[[i]]]
  }
  
}, error = function(e){print(paste0("SGL cannot fit"))})
fit_table[, "SGL"] <- result_compare(Y_fit_SGL[observed_index], Y[observed_index])[[1]]

save(list = c('results_LOO', 
              'results_LOO_MSE',
              "fit_table",
              "results_feature_selection_list"), file = "Results_LOO_FS.RData")