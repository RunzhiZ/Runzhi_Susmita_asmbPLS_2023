######## PARAMETERS ########**********************************************
#### n: number of samples
#### c: censoring rate
#### r: noise to signal ratio
#### predictor: predictor matrix
#### beta: beta for predictor matrix
#### seed: seed for data generation
####**********************************************************************

#### Different strategies for beta ####
## First 5 - 1:5,  remaining 0
f <- 5
beta_1_microbiome <- c(1:5, rep(0, q - f))
beta_1_metabolome <- c(1:5, rep(0, p - f))
beta_1 <- c(beta_1_microbiome, beta_1_metabolome)
## beta = exp(-j)
beta_2_microbiome <- exp(- (1:q))
beta_2_metabolome <- exp(- (1:p))
beta_2 <- c(beta_2_microbiome, beta_2_metabolome)
## beta = 1/j
beta_3_microbiome <- 1 / (1:q)
beta_3_metabolome <- 1 / (1:p)
beta_3 <- c(beta_3_microbiome, beta_3_metabolome)
## beta = 1 for all
beta_4_microbiome <- rep(1, q)
beta_4_metabolome <- rep(1, p)
beta_4 <- c(beta_4_microbiome, beta_4_metabolome)

#### Survival time ####
survival_time_simulation <- function(n, n_censored, c, predictor, beta, r, seed) {
  lognormal_survival <- weibull_survival <- matrix(NA, nrow = n, ncol = 4)
  colnames(lognormal_survival) <- colnames(weibull_survival) <- c("log_Survival_time", "log_Censored_time","log_right_censored_time", "Event_indicator")
  beta_matrix <- matrix(beta, ncol = 1)
  cov_predictor_train <- cov(predictor[1:n_censored,])
  cov_predictor_test <- cov(predictor[(n_censored+1):n,])

  sigma_sqr_train <- as.vector(t(beta_matrix) %*% cov_predictor_train %*% beta_matrix)
  sigma_sqr_test <- as.vector(t(beta_matrix) %*% cov_predictor_test %*% beta_matrix)
  
  set.seed(seed)
  ## T - lognormal
  error_lognormal <- c(sqrt(r) * sqrt(sigma_sqr_train) * rnorm(n_censored), sqrt(r) * sqrt(sigma_sqr_test) * rnorm(n-n_censored))
  logT_lognormal <- predictor %*% beta_matrix + error_lognormal + log(500)
  lognormal_survival[, 1] <- logT_lognormal
  ## T - Weibull
  weibull_5_1_train <- rweibull(n_censored, 5, 1)
  weibull_5_1_test <- rweibull(n - n_censored, 5, 1)
  error_weibull <- c(sqrt(r) * sqrt(sigma_sqr_train) * (log(weibull_5_1_train) - mean(log(weibull_5_1_train))) / sqrt(var(log(weibull_5_1_train))), 
                     sqrt(r) * sqrt(sigma_sqr_test) * (log(weibull_5_1_test) - mean(log(weibull_5_1_test))) / sqrt(var(log(weibull_5_1_test))))
  logT_weibull <- predictor %*% beta_matrix + error_weibull + log(500)
  weibull_survival[, 1] <- logT_weibull
  ## censoring rate
  rate <- qnorm(c,lower.tail = F)
  cencoring <- rnorm(n_censored, rate * sqrt(sigma_sqr_train * (1 + r)), sqrt(sigma_sqr_train * (1 + r)))
  lognormal_survival[, 2] <- logT_lognormal + c(cencoring, rep(1,n - n_censored))
  weibull_survival[, 2] <- logT_weibull + c(cencoring, rep(1,n - n_censored))
  lognormal_survival[, 3] <- apply(lognormal_survival[,1:2],1,min)
  weibull_survival[, 3] <- apply(weibull_survival[,1:2],1,min)
  lognormal_survival[, 4] <- as.numeric(lognormal_survival[, 1] <= lognormal_survival[, 2])
  weibull_survival[, 4] <- as.numeric(weibull_survival[, 1] <= weibull_survival[, 2])
  return(list(lognormal_survival = lognormal_survival, 
              weibull_survival = weibull_survival, 
              sigma_sqr_test = sigma_sqr_test))
}

