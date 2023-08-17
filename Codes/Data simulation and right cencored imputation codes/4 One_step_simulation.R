######## PARAMETERS ########**********************************************
#### n: number of samples
#### q: number of microbial features
#### p: number of metabolome features
#### m: number of baseline sequence reads for microbial data
#### n_effect: number of pairs of effect between microbial and metabolome features
#### c: censoring rate
#### r: noise to signal ratio
#### beta: beta for predictor matrix
#### seed: seed for data generation
####**********************************************************************

#### scale function, for column with all 0, keep the original column ####
X_scale_0 <- function(data) {
  if (sum(data)!=0) {
    data <- scale(data)
  }
  return(data)
}
#### compare function: compare predicted and true Y ####
result_compare <- function (prediction, true) {
  result_abs = sum(abs(prediction-true))/length(true)
  result_sq = sum((prediction-true)^2)/length(true)
  return(list(result_sq = result_sq,
              result_abs = result_abs))
}
#### ONE STEP SIMULATION FUNCTION ####
one_step_simulation <- function(n, n_censored, q, p, m, n_effect, effect_scale, r, c, beta, seed, pseudo_impute = T) {
  microbiome_count <- microbiome_data_generate(n, q, m, seed)$microbiome_raw_data
  
  microbiome_relative_abundance <- microbiome_data_generate(n, q, m, seed)$microbiome_relative_abundance
  microbiome_scale <- apply(microbiome_relative_abundance, 2, X_scale_0)
  metabolome_data <- metabolome_data_generate(n, n_effect, q, p, effect_scale, microbiome_scale, seed = seed)$metabolome_data
  metabolome_scale <- apply(metabolome_data, 2, X_scale_0)
  
  if(pseudo_impute == T){
    microbiome_count[microbiome_count == 0] <- 0.5
    microbiome_relative_abundance_pseudo <- t(apply(microbiome_count, 1, function(data) {
      data / sum(data)
    }))
    predictor_data <- cbind(microbiome_relative_abundance_pseudo, metabolome_data)
  } else {
    predictor_data <- cbind(microbiome_relative_abundance, metabolome_data)
  }
  
  predictor_data_scale <- cbind(microbiome_scale, metabolome_scale)
  survival_time_lognormal <- survival_time_simulation(n, n_censored, c, predictor_data_scale, beta, r, seed)$lognormal_survival
  survival_time_weibull <- survival_time_simulation(n, n_censored, c, predictor_data_scale, beta, r, seed)$weibull_survival
  sigma_sqr <- survival_time_simulation(n, n_censored, c, predictor_data_scale, beta, r, seed)$sigma_sqr_test
  return(list(predictor_data = predictor_data,
              survival_time_lognormal = survival_time_lognormal,
              survival_time_weibull = survival_time_weibull,
              sigma_sqr = sigma_sqr
  ))
}

## example ##
n = 200
n_censored = 100
q = 200
p = 200
m = 20000
n_effect = 100
effect_scale = 0.5
r = 1
c = 0.3
beta = beta_1
seed = 123

one_step_simulation(n = 200, n_censored = 100, q = 200, p = 200, m = 2000,
                    n_effect = 100, effect_scale = 0.5,
                    r = 0.1, c = 0.3, beta = beta_1,
                    seed = 1)
