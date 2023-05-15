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
library(asmbPLS)

setwd("/blue/datta/runzhi.zhang/asmbPLS/")
load("Cencored data imputation function.RData")
load("Data Simulation function.RData")

###############################################################################
# parameters setting
###############################################################################
dist <- "lognormal" ## or "weibull"
imputation_method <- "mean_imputation"

## fixed part ##
iter = 100
n = 200
n_censored = 100
k = 5
PLS_term_selected = 5

## different options
# low dimension
q = 200
p = 200
m = 20000
n_effect = 200

# mixed dimension
q = 1000
p = 200
m = 100000
n_effect = 500

# high dimension
q = 1000
p = 1000
m = 100000
n_effect = 1200

## tuning ##
c = 0 ## or c = 0.1, 0.3, 0.5, 0.6, 0.7
r = 0 ## or r = 0.1, 0.2, 0.5, 1
effect_scale = 0.5 # or = 2

## dimension for blocks ##
X_dim <- c(q, p)

#### lambda combination ####
lambda_percent <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.975)
lambda_percent_0.5 <- c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.975)
lambda_percent_0.7 <- c(0.7, 0.8, 0.9, 0.95, 0.975)
lambda_percent_table <- cbind(rep(lambda_percent, each = length(lambda_percent)), rep(lambda_percent, length(lambda_percent)))
lambda_percent_table_0.5 <- cbind(rep(lambda_percent_0.5, each = length(lambda_percent_0.5)), rep(lambda_percent_0.5, length(lambda_percent_0.5)))
lambda_percent_table_0.7 <- cbind(rep(lambda_percent_0.7, each = length(lambda_percent_0.7)), rep(lambda_percent_0.7, length(lambda_percent_0.7)))

#### Different strategies for beta ####
## First 5 - ,1:5,  remaining 0
f <- 5
beta_1_microbiome <- c(1:5, rep(0, q - f))
beta_1_metabolome <- c(1:5, rep(0, p - f))
beta_1 <- c(beta_1_microbiome, beta_1_metabolome)
beta_1 = beta_1/c(sqrt(t(beta_1)%*%beta_1))
## beta = exp(-j)
beta_2_microbiome <- exp(- (1:q))
beta_2_metabolome <- exp(- (1:p))
beta_2 <- c(beta_2_microbiome, beta_2_metabolome)
beta_2 = beta_2/c(sqrt(t(beta_2)%*%beta_2))
## beta = 1/j
beta_3_microbiome <- 1 / (1:q)
beta_3_metabolome <- 1 / (1:p)
beta_3 <- c(beta_3_microbiome, beta_3_metabolome)
beta_3 = beta_3/c(sqrt(t(beta_3)%*%beta_3))
## beta = 1 for all
beta_4_microbiome <- rep(1, q)
beta_4_metabolome <- rep(1, p)
beta_4 <- c(beta_4_microbiome, beta_4_metabolome)
beta_4 = beta_4/c(sqrt(t(beta_4)%*%beta_4))
## different proportion for different blocks
beta_5_microbiome <- c(1:5, 1:5, rep(0, q - 10))
beta_5_metabolome <- c(1:5, rep(0, p - 5))
beta_5 <- c(beta_5_microbiome, beta_5_metabolome)
beta_5 = beta_5/c(sqrt(t(beta_5)%*%beta_5))
## different proportion for different blocks
beta_6_microbiome <- c(1:5, rep(0, q - 5))
beta_6_metabolome <- c(1:5, 1:5, rep(0, p - 10))
beta_6 <- c(beta_6_microbiome, beta_6_metabolome)
beta_6 = beta_6/c(sqrt(t(beta_6)%*%beta_6))

###############################################################################
# Tables to save all the results:
# 1) MSE of different methods
# 2) computation time of different methods
###############################################################################
#### result table ####
results_feature_selection_list <- list()
results_asmbPLS_predict_list <- results_asmbPLS_fit_list <- results_asmbPLS_cv_list <- results_asmbPLS_predict_outlier_list <- list()
results_asmbPLS_0.5_predict_list <- results_asmbPLS_0.5_fit_list <- results_asmbPLS_0.5_cv_list <- results_asmbPLS_0.5_predict_outlier_list <- list()
results_asmbPLS_0.7_predict_list <- results_asmbPLS_0.7_fit_list <- results_asmbPLS_0.7_cv_list <- results_asmbPLS_0.7_predict_outlier_list <- list()
results_mbPLS_predict_list <- results_mbPLS_fit_list <- list()
results_BF_predict_list <- results_BF_fit_list <- list()
results_IPF_predict_list <- results_IPF_fit_list <- list()
results_SGL_predict_list <- results_SGL_fit_list <- list()
results_priority_lasso_predict_list <- results_priority_lasso_fit_list <- list()

results_simulation_sigma_sqr_list <- list()

#### time table ####
time_table_1 <- 
  time_table_2 <- 
  time_table_3 <- 
  time_table_4 <- 
  time_table_5 <- 
  time_table_6 <- matrix(NA, nrow = iter, ncol = 9)
colnames(time_table_1) <- 
  colnames(time_table_2) <- 
  colnames(time_table_3) <- 
  colnames(time_table_4) <- 
  colnames(time_table_5) <- 
  colnames(time_table_6) <- 
  c("imputation", "asmbPLS", "asmbPLS_0.5", "asmbPLS_0.7", "mbPLS", 
    "BF", "IPF", "SGL", "priority_lasso")
time_table_1[, "imputation"] <- 
  time_table_2[, "imputation"] <- 
  time_table_3[, "imputation"] <- 
  time_table_4[, "imputation"] <- 
  time_table_5[, "imputation"] <- 
  time_table_6[, "imputation"] <- imputation_method

###############################################################################
# Simulation and comparison:
# 1) asmbPLS
# 2) asmbPLS_0.5
# 3) asmbPLS_0.7
# 4) mbPLS (only prediction, without feature selection)
# 5) block forest
# 6) IPF-lasso
# 7) SGL
# 8) priority-lasso
###############################################################################

for (b in 1:6) {
  print(paste0("****************************** Beta ", b, " ******************************"))
  #### generate matrix to save the results ####
  eval(parse(text = paste0("results_feature_selection_list[[", b, "]] <- list()")))
  
  eval(parse(text = paste0("results_", "asmbPLS", "_predict_list[[", b, "]] <- matrix(NA, nrow = PLS_term_selected, ncol = iter)")))
  eval(parse(text = paste0("results_", "asmbPLS", "_cv_list[[", b, "]] <- matrix(NA, nrow = PLS_term_selected, ncol = iter)")))
  eval(parse(text = paste0("results_", "asmbPLS", "_fit_list[[", b, "]] <- matrix(NA, nrow = PLS_term_selected, ncol = iter)")))
  eval(parse(text = paste0("results_", "asmbPLS", "_predict_outlier_list[[", b, "]] <- matrix(NA, nrow = PLS_term_selected, ncol = iter)")))
  
  eval(parse(text = paste0("results_", "asmbPLS", "_0.5_predict_list[[", b, "]] <- matrix(NA, nrow = PLS_term_selected, ncol = iter)")))
  eval(parse(text = paste0("results_", "asmbPLS", "_0.5_cv_list[[", b, "]] <- matrix(NA, nrow = PLS_term_selected, ncol = iter)")))
  eval(parse(text = paste0("results_", "asmbPLS", "_0.5_fit_list[[", b, "]] <- matrix(NA, nrow = PLS_term_selected, ncol = iter)")))
  eval(parse(text = paste0("results_", "asmbPLS", "_0.5_predict_outlier_list[[", b, "]] <- matrix(NA, nrow = PLS_term_selected, ncol = iter)")))
  
  eval(parse(text = paste0("results_", "asmbPLS", "_0.7_predict_list[[", b, "]] <- matrix(NA, nrow = PLS_term_selected, ncol = iter)")))
  eval(parse(text = paste0("results_", "asmbPLS", "_0.7_cv_list[[", b, "]] <- matrix(NA, nrow = PLS_term_selected, ncol = iter)")))
  eval(parse(text = paste0("results_", "asmbPLS", "_0.7_fit_list[[", b, "]] <- matrix(NA, nrow = PLS_term_selected, ncol = iter)")))
  eval(parse(text = paste0("results_", "asmbPLS", "_0.7_predict_outlier_list[[", b, "]] <- matrix(NA, nrow = PLS_term_selected, ncol = iter)")))
  
  eval(parse(text = paste0("results_", "mbPLS", "_predict_list[[", b, "]] <- matrix(NA, nrow = PLS_term_selected, ncol = iter)")))
  eval(parse(text = paste0("results_", "mbPLS", "_fit_list[[", b, "]] <- matrix(NA, nrow = PLS_term_selected, ncol = iter)")))
  
  eval(parse(text = paste0("results_", "BF", "_predict_list[[", b, "]] <- numeric()")))
  eval(parse(text = paste0("results_", "BF", "_fit_list[[", b, "]] <- numeric()")))
  
  eval(parse(text = paste0("results_", "IPF", "_predict_list[[", b, "]] <- numeric()")))
  eval(parse(text = paste0("results_", "IPF", "_fit_list[[", b, "]] <- numeric()")))
  
  eval(parse(text = paste0("results_", "SGL", "_predict_list[[", b, "]] <- numeric()")))
  eval(parse(text = paste0("results_", "SGL", "_fit_list[[", b, "]] <- numeric()")))
  
  eval(parse(text = paste0("results_", "priority_lasso", "_predict_list[[", b, "]] <- numeric()")))
  eval(parse(text = paste0("results_", "priority_lasso", "_fit_list[[", b, "]] <- numeric()")))
  
  eval(parse(text = paste0("results_", "simulation", "_sigma_sqr_list[[", b, "]] <- numeric()")))
  
  beta <- eval(parse(text = paste0("beta_", b)))
  
  for (l in 1:iter) {
    print(paste0("************************** Iteration ", l, " **************************"))
    survival_data <- one_step_simulation(n = n, n_censored = n_censored, q = q, p = p, m = m,
                                         n_effect = n_effect, effect_scale = effect_scale,
                                         r = r, c = c, beta = beta,
                                         seed = l)
    eval(parse(text = paste0("results_", "simulation", "_sigma_sqr_list[[", b, "]][l] <- survival_data$sigma_sqr")))
    
    
    #### data pre-processing ####
    X <- survival_data$predictor_data
    ## imputation method
    if (imputation_method == "mean_imputation") {
      if (c > 0) {
        Y_train <- matrix((mean_imputation(n_censored, survival_data$survival_time_lognormal[1:n_censored, 3:4], F)$imputed_table[, 4]))
      } else {
        Y_train <- matrix((survival_data$survival_time_lognormal[1:n_censored, 3]))
      }
    }
    if (imputation_method == "reweighting") {
      if (c > 0) {
        Y_train <- matrix((reweighting(n_censored, survival_data$survival_time_lognormal[1:n_censored, 3:4], F)$imputed_table[, 4]))
      } else {
        Y_train <- matrix((survival_data$survival_time_lognormal[1:n_censored, 3]))
      }
    }
    ### clr transformation for X ####
    X_clr <- X
    X_clr[, 1:X_dim[1]] <- clr(X[, 1:X_dim[1]])
    train_index <- 1:n_censored
    X_train <- X_clr[train_index, ]
    X_test <- X_clr[-train_index, ]
    Y_test <- as.matrix((survival_data$survival_time_lognormal[(n_censored+1):n, 3]))
    Y_indicator <- survival_data$survival_time_lognormal[1:n_censored, 4]
    
    #### blocks setting ####
    blocks_vector <- rep(1:length(X_dim), times = X_dim)
    blocks <- lapply(1:2, function(x) which(blocks_vector == x))
    blocks_1 <- lapply(1:2, function(x) which(blocks_vector == x))
    blocks_2 <- lapply(2:1, function(x) which(blocks_vector == x))
    
    #### features significance ####
    results_feature_selection_list[[b]][[l]] <- list()
    for(i in 1:length(X_dim)) {
      results_feature_selection_list[[b]][[l]][[i]] <- matrix(NA, nrow = X_dim[i], ncol = 10)
      row.names(results_feature_selection_list[[b]][[l]][[i]]) <- colnames(X_train)[blocks[[i]]]
      colnames(results_feature_selection_list[[b]][[l]][[i]]) <- c("p_value", "p_value_adjusted", "asmbPLS", "asmbPLS_0.5", "asmbPLS_0.7", "mbPLS", "BF", "IPF", "SGL", "priority_lasso")
    }
    for(i in 1:length(X_dim)) {
      for(j in 1:X_dim[i]) {
        if(sum(X_train[,blocks[[i]][j]]!=0)) {results_feature_selection_list[[b]][[l]][[i]][j, "p_value"] <- summary(lm(Y_train ~ X_train[,blocks[[i]][j]]))$coefficients[2, 4]}
      }
      results_feature_selection_list[[b]][[l]][[i]][, "p_value_adjusted"] <- p.adjust(results_feature_selection_list[[b]][[l]][[i]][, "p_value"], method = "BH")
    }
    
    
    #### asmbPLS ####
    time_start <- Sys.time()
    asmbPLS_results_cv <- asmbPLS.cv(X.matrix = X_train, 
                                     Y.matrix = Y_train,
                                     PLS.comp = PLS_term_selected, 
                                     X.dim = X_dim, 
                                     quantile.comb.table = lambda_percent_table, 
                                     Y.indicator = Y_indicator,
                                     k = k, 
                                     ncv = 1)
    PLS_table <- asmbPLS_results_cv$quantile_table_CV[,1:length(X_dim)]
    asmbPLS_results <- asmbPLS.fit(X.matrix = X_train, 
                                   Y.matrix = Y_train, 
                                   PLS.comp = PLS_term_selected, 
                                   X.dim = X_dim,
                                   quantile.comb = PLS_table[1:PLS_term_selected,])
    time_end <- Sys.time()
    time_asmbPLS <- difftime(time_end, time_start, units = "secs")
    print(paste0("Time for asmbPLS model fit: ", round(time_asmbPLS,4)))
    results_asmbPLS_predict_single <- NULL
    results_asmbPLS_fit_single <- NULL
    for (PLS_term in 1:PLS_term_selected) {
      Y_predict_asmbPLS <- asmbPLS.predict(fit.results = asmbPLS_results,
                                           X.matrix.new = X_test, 
                                           PLS.comp = PLS_term)
      Y_fit_asmbPLS <- asmbPLS.predict(fit.results = asmbPLS_results,
                                       X.matrix.new = X_train, 
                                       PLS.comp = PLS_term)
      results_asmbPLS_predict_single[PLS_term] <- result_compare(Y_predict_asmbPLS$Y_pred, Y_test)[[1]]
      results_asmbPLS_fit_single[PLS_term] <- result_compare(Y_fit_asmbPLS$Y_pred, Y_train)[[1]]
    }
    eval(parse(text = paste0("results_", "asmbPLS", "_predict_list[[", b, "]][, l] <- results_asmbPLS_predict_single")))
    eval(parse(text = paste0("results_", "asmbPLS", "_cv_list[[", b, "]][, l] <- asmbPLS_results_cv$min_result_table[, 3]")))
    eval(parse(text = paste0("results_", "asmbPLS", "_fit_list[[", b, "]][, l] <- results_asmbPLS_fit_single")))
    
    ## feature selection
    for(i in 1:length(X_dim)) {
      results_feature_selection_list[[b]][[l]][[i]][, "asmbPLS"] <- asmbPLS_results$X_weight[[i]][,1]
    }
    
    
    #### asmbPLS_0.5 ####
    time_start <- Sys.time()
    asmbPLS_results_cv <- asmbPLS.cv(X.matrix = X_train, 
                                     Y.matrix = Y_train,
                                     PLS.comp = PLS_term_selected, 
                                     X.dim = X_dim, 
                                     quantile.comb.table = lambda_percent_table_0.5, 
                                     Y.indicator = Y_indicator,
                                     k = k, 
                                     ncv = 1)
    PLS_table <- asmbPLS_results_cv$quantile_table_CV[,1:length(X_dim)]
    asmbPLS_results <- asmbPLS.fit(X.matrix = X_train, 
                                   Y.matrix = Y_train, 
                                   PLS.comp = PLS_term_selected, 
                                   X.dim = X_dim,
                                   quantile.comb = PLS_table[1:PLS_term_selected,])
    time_end <- Sys.time()
    time_asmbPLS_0.5 <- difftime(time_end, time_start, units = "secs")
    print(paste0("Time for asmbPLS_0.5 model fit: ", round(time_asmbPLS_0.5,4)))
    results_asmbPLS_predict_single <- NULL
    results_asmbPLS_fit_single <- NULL
    for (PLS_term in 1:PLS_term_selected) {
      Y_predict_asmbPLS <- asmbPLS.predict(fit.results = asmbPLS_results,
                                           X.matrix.new = X_test, 
                                           PLS.comp = PLS_term)
      Y_fit_asmbPLS <- asmbPLS.predict(fit.results = asmbPLS_results,
                                       X.matrix.new = X_train, 
                                       PLS.comp = PLS_term)
      results_asmbPLS_predict_single[PLS_term] <- result_compare(Y_predict_asmbPLS$Y_pred, Y_test)[[1]]
      results_asmbPLS_fit_single[PLS_term] <- result_compare(Y_fit_asmbPLS$Y_pred, Y_train)[[1]]
    }
    eval(parse(text = paste0("results_", "asmbPLS", "_0.5_predict_list[[", b, "]][, l] <- results_asmbPLS_predict_single")))
    eval(parse(text = paste0("results_", "asmbPLS", "_0.5_cv_list[[", b, "]][, l] <- asmbPLS_results_cv$min_result_table[, 3]")))
    eval(parse(text = paste0("results_", "asmbPLS", "_0.5_fit_list[[", b, "]][, l] <- results_asmbPLS_fit_single")))
    
    ## feature selection
    for(i in 1:length(X_dim)) {
      results_feature_selection_list[[b]][[l]][[i]][, "asmbPLS_0.5"] <- asmbPLS_results$X_weight[[i]][,1]
    }
    
    
    #### asmbPLS_0.7 ####
    time_start <- Sys.time()
    asmbPLS_results_cv <- asmbPLS.cv(X.matrix = X_train, 
                                     Y.matrix = Y_train,
                                     PLS.comp = PLS_term_selected, 
                                     X.dim = X_dim, 
                                     quantile.comb.table = lambda_percent_table_0.7, 
                                     Y.indicator = Y_indicator,
                                     k = k, 
                                     ncv = 1)
    PLS_table <- asmbPLS_results_cv$quantile_table_CV[,1:length(X_dim)]
    asmbPLS_results <- asmbPLS.fit(X.matrix = X_train, 
                                   Y.matrix = Y_train, 
                                   PLS.comp = PLS_term_selected, 
                                   X.dim = X_dim,
                                   quantile.comb = PLS_table[1:PLS_term_selected,])
    time_end <- Sys.time()
    time_asmbPLS_0.7 <- difftime(time_end, time_start, units = "secs")
    print(paste0("Time for asmbPLS_0.7 model fit: ", round(time_asmbPLS_0.7,4)))
    results_asmbPLS_predict_single <- NULL
    results_asmbPLS_fit_single <- NULL
    for (PLS_term in 1:PLS_term_selected) {
      Y_predict_asmbPLS <- asmbPLS.predict(fit.results = asmbPLS_results,
                                           X.matrix.new = X_test, 
                                           PLS.comp = PLS_term)
      Y_fit_asmbPLS <- asmbPLS.predict(fit.results = asmbPLS_results,
                                       X.matrix.new = X_train, 
                                       PLS.comp = PLS_term)
      results_asmbPLS_predict_single[PLS_term] <- result_compare(Y_predict_asmbPLS$Y_pred, Y_test)[[1]]
      results_asmbPLS_fit_single[PLS_term] <- result_compare(Y_fit_asmbPLS$Y_pred, Y_train)[[1]]
    }
    eval(parse(text = paste0("results_", "asmbPLS", "_0.7_predict_list[[", b, "]][, l] <- results_asmbPLS_predict_single")))
    eval(parse(text = paste0("results_", "asmbPLS", "_0.7_cv_list[[", b, "]][, l] <- asmbPLS_results_cv$min_result_table[, 3]")))
    eval(parse(text = paste0("results_", "asmbPLS", "_0.7_fit_list[[", b, "]][, l] <- results_asmbPLS_fit_single")))
    
    ## feature selection
    for(i in 1:length(X_dim)) {
      results_feature_selection_list[[b]][[l]][[i]][, "asmbPLS_0.7"] <- asmbPLS_results$X_weight[[i]][,1]
    }
    
    
    #### mbPLS ####
    time_start <- Sys.time()
    mbPLS_results <- mbPLS_fit_rcpp(X.matrix = X_train, 
                                    Y.matrix = Y_train, 
                                    PLS.comp = PLS_term_selected, 
                                    X.dim = X_dim)
    time_end <- Sys.time()
    time_mbPLS <- difftime(time_end, time_start, units = "secs")
    print(paste0("Time for mbPLS model fit: ", round(time_mbPLS,4)))
    results_mbPLS_predict_single <- NULL
    results_mbPLS_fit_single <- NULL
    for (PLS_term in 1:PLS_term_selected) {
      Y_predict_mbPLS <- asmbPLS.predict(fit.results = mbPLS_results,
                                           X.matrix.new = X_test, 
                                           PLS.comp = PLS_term)
      Y_fit_mbPLS <- asmbPLS.predict(fit.results = mbPLS_results,
                                       X.matrix.new = X_train, 
                                       PLS.comp = PLS_term)
      results_mbPLS_predict_single[PLS_term] <- result_compare(Y_predict_mbPLS$Y_pred, Y_test)[[1]]
      results_mbPLS_fit_single[PLS_term] <- result_compare(Y_fit_mbPLS$Y_pred, Y_train)[[1]]
    }
    eval(parse(text = paste0("results_", "mbPLS", "_predict_list[[", b, "]][, l] <- results_mbPLS_predict_single")))
    eval(parse(text = paste0("results_", "mbPLS", "_fit_list[[", b, "]][, l] <- results_mbPLS_fit_single")))
    
    ## feature selection
    for(i in 1:length(X_dim)) {
      results_feature_selection_list[[b]][[l]][[i]][, "mbPLS"] <- mbPLS_results$x_weight[[i]][,1]
    }
    
    
    #### Block forest ####
    Y_train_vector <- as.vector(Y_train)
    time_start <- Sys.time()
    blockforobj <- blockfor(X_train, Y_train_vector, replace = TRUE, blocks = blocks, block.method = "BlockForest", importance = "impurity", seed = l)
    time_end <- Sys.time()
    time_BF <- difftime(time_end, time_start, units = "secs")
    print(paste0("Time for BF model fit: ", round(time_BF,4)))
    Y_predict_BF <- predict(blockforobj$forest, data = X_test)
    Y_fit_BF <- predict(blockforobj$forest, data = X_train)
    eval(parse(text = paste0("results_", "BF", "_predict_list[[", b, "]][l] <- result_compare(Y_predict_BF$predictions, Y_test)[[1]]")))
    eval(parse(text = paste0("results_", "BF", "_fit_list[[", b, "]][l] <- result_compare(Y_fit_BF$predictions, Y_train)[[1]]")))
    
    ## feature selection
    for(i in 1:length(X_dim)) {
      results_feature_selection_list[[b]][[l]][[i]][, "BF"] <- blockforobj$forest$variable.importance[blocks[[i]]]
    }
    
    
    #### IPF ####
    set.seed(l)
    Y_train_vector <- as.vector(Y_train)
    pflist <- list(c(1, 1), c(2, 1), c(1, 2), c(3, 1), c(1, 3), c(4, 1), c(1, 4), c(5, 1), c(1, 5))
    time_start <- Sys.time()
    ipf_fit <- cvr2.ipflasso(X = X_train, Y = Y_train_vector, family = "gaussian", type.measure = "mse",
                             blocks = blocks,
                             pflist = pflist, nfolds = 5, ncv = 10)
    time_end <- Sys.time()
    time_IPF <- difftime(time_end, time_start, units = "secs")
    print(paste0("Time for IPF model fit: ", round(time_IPF,4)))
    Y_predict_IPF <- ipflasso.predict(ipf_fit, X_test)
    Y_fit_IPF <- ipflasso.predict(ipf_fit, X_train)
    eval(parse(text = paste0("results_", "IPF", "_predict_list[[", b, "]][l] <- result_compare(Y_predict_IPF$linpredtest, Y_test)[[1]]")))
    eval(parse(text = paste0("results_", "IPF", "_fit_list[[", b, "]][l] <- result_compare(Y_fit_IPF$linpredtest, Y_train)[[1]]")))
    
    ## feature selection
    for(i in 1:length(X_dim)) {
      results_feature_selection_list[[b]][[l]][[i]][, "IPF"] <- ipf_fit$coeff[,ipf_fit$ind.bestlambda][-1][blocks[[i]]]
    }
    
    
    #### SGL ####
    tryCatch({
      data = list(x = X_train, y = Y_train)
      index <- c(rep(1, X_dim[1]), rep(2, X_dim[2]))
      time_start <- Sys.time()
      cvFit = cvSGL(data, index, type = "linear", nfold = 5)
      Fit = SGL(data, index, type = "linear", lambdas = cvFit$lambdas)
      time_end <- Sys.time()
      time_SGL <- difftime(time_end, time_start, units = "secs")
      print(paste0("Time for SGL model fit: ", round(time_SGL,4)))
      Y_predict_SGL <- predictSGL(Fit, X_test, which.min(cvFit$lldiff))
      Y_fit_SGL <- predictSGL(Fit, X_train, which.min(cvFit$lldiff))
      eval(parse(text = paste0("results_", "SGL", "_predict_list[[", b, "]][l] <- result_compare(Y_predict_SGL, Y_test)[[1]]")))
      eval(parse(text = paste0("results_", "SGL", "_fit_list[[", b, "]][l] <- result_compare(Y_fit_SGL, Y_train)[[1]]")))
      
      ## feature selection
      for(i in 1:length(X_dim)) {
        results_feature_selection_list[[b]][[l]][[i]][, "SGL"] <- Fit$beta[, which.min(cvFit$lldiff)][blocks[[i]]]
      }
      
      eval(parse(text = paste0("time_table_", b, "[l, 8] <- round(time_SGL, 4)")))
      
    }, error = function(e){print(paste0("SGL, beta_", b, ": itertation ",l, " cannot fit"))})
    
    
    #### priority lasso ####
    time_start <- Sys.time()
    priority_fit <- cvm_prioritylasso(X = X_train, Y = Y_train, family = "gaussian", 
                                      type.measure = "mse", standardize = T, 
                                      blocks.list = list(blocks_1, blocks_2), block1.penalization = TRUE, lambda.type = "lambda.1se")
    
    time_end <- Sys.time()
    time_priority_lasso <- difftime(time_end, time_start, units = "secs")
    print(paste0("Time for priority lasso model fit: ", round(time_priority_lasso,4)))
    Y_predict_priority_lasso <- predict(priority_fit, X_test, type = "response")
    Y_fit_priority_lasso <- predict(priority_fit, X_train, type = "response")
    eval(parse(text = paste0("results_", "priority_lasso", "_predict_list[[", b, "]][l] <- result_compare(Y_predict_priority_lasso, Y_test)[[1]]")))
    eval(parse(text = paste0("results_", "priority_lasso", "_fit_list[[", b, "]][l] <- result_compare(Y_fit_priority_lasso, Y_train)[[1]]")))
    
    ## feature selection
    for(i in 1:length(X_dim)) {
      results_feature_selection_list[[b]][[l]][[i]][, "priority_lasso"] <- priority_fit$coefficients[blocks[[i]]]
    }
    
    #### time ####
    eval(parse(text = paste0("time_table_", b, "[l, 2] <- round(time_asmbPLS, 4)")))
    eval(parse(text = paste0("time_table_", b, "[l, 3] <- round(time_asmbPLS_0.5, 4)")))
    eval(parse(text = paste0("time_table_", b, "[l, 4] <- round(time_asmbPLS_0.7, 4)")))
    eval(parse(text = paste0("time_table_", b, "[l, 5] <- round(time_mbPLS, 4)")))
    eval(parse(text = paste0("time_table_", b, "[l, 6] <- round(time_BF, 4)")))
    eval(parse(text = paste0("time_table_", b, "[l, 7] <- round(time_IPF, 4)")))
    eval(parse(text = paste0("time_table_", b, "[l, 9] <- round(time_priority_lasso, 4)")))
  }
}

###############################################################################
# Save results
###############################################################################
path <- paste0("/blue/datta/runzhi.zhang/asmbPLS/q.", q, ".p.", p, "/", dist, "_", imputation_method, "/Results_V4/")
## r
if(r == 0) {r_index <- "00"}
if(r == 0.1) {r_index <- "01"}
if(r == 0.2) {r_index <- "02"}
if(r == 0.5) {r_index <- "05"}
if(r == 1) {r_index <- "10"}
## c
c_index <- c*10
## effect_scale
if(effect_scale == 0.5) {e_index <- "05"}
if(effect_scale == 2) {e_index <- "20"}

file_name <- paste0("Results_c0", c_index, "_r", r_index, "_e", e_index, "_", dist, ".RData")
setwd(path)
save(list = c('results_asmbPLS_predict_list',
              'results_asmbPLS_predict_outlier_list',
              "results_asmbPLS_cv_list",
              'results_asmbPLS_fit_list',
              
              'results_asmbPLS_0.5_predict_list',
              'results_asmbPLS_0.5_predict_outlier_list',
              "results_asmbPLS_0.5_cv_list",
              'results_asmbPLS_0.5_fit_list',
              
              'results_asmbPLS_0.7_predict_list',
              'results_asmbPLS_0.7_predict_outlier_list',
              "results_asmbPLS_0.7_cv_list",
              'results_asmbPLS_0.7_fit_list',
              
              'results_mbPLS_predict_list',
              'results_mbPLS_fit_list',
              
              'results_BF_predict_list',
              'results_BF_fit_list',
              
              'results_IPF_predict_list',
              'results_IPF_fit_list',
              
              'results_SGL_predict_list',
              'results_SGL_fit_list',
              
              'results_priority_lasso_predict_list',
              'results_priority_lasso_fit_list',
              
              
              "results_simulation_sigma_sqr_list",
              
              "results_feature_selection_list",
              
              "time_table_1",
              "time_table_2",
              "time_table_3",
              "time_table_4",
              "time_table_5",
              "time_table_6"
),
file = file_name)
