######## PARAMETERS ########**********************************************
#### n: number of samples
#### q: number of microbial features
#### p: number of metabolome features
#### n_effect: number of pairs of effect between microbial and metabolome features
#### seed: seed for data generation
####**********************************************************************

#### Microbial effect on metabolome ####
effect_index <- function(n_effect, q, p, seed) {
  comb <- p*q
  index_table <- cbind(rep(1:q, p), rep(1:p, each = q))
  set.seed(seed)
  index_table_return <- index_table[sample(1:comb, n_effect, replace = F), ]
  return(index_table_return)
}

metabolome_data_generate <- function(n, n_effect, q, p, effect_scale, microbiome_matrix, base_scale = 4, error_mean = 0, error_sd = 1, seed) {
  ## Microbial effect
  set.seed(seed)
  effect_matrix <- matrix(0, nrow = q, ncol = p)
  index_table <- effect_index(n_effect, q, p, seed)
  for (i in 1:nrow(index_table)) {
    ## half positive, half negative
    if (i %% 2 == 1) {
      effect <- runif(1, effect_scale, 2 * effect_scale)
    } else {
      effect <- runif(1, -2 * effect_scale, -1 * effect_scale)
    }
    effect_matrix[index_table[i, 1], index_table[i, 2]] <- effect
  }
  micro_effect <- microbiome_matrix %*% effect_matrix
  ## Metabolome base
  metabolome_base <- metabolome_error <- matrix(NA, nrow = n, ncol = p)
  for (i in 1:ncol(metabolome_base)) {
    metabolome_base[, i] <- runif(1, base_scale, 2 * base_scale)
  }
  ## Metabolome error
  for (i in 1:p) {
    metabolome_error[, i] <- rnorm(n, error_mean, error_sd)
  }
  ## Metabolome data
  metabolome_data <- metabolome_base + micro_effect + metabolome_error
  metabolome_data[(metabolome_data <= 0)] <- 0
  colnames(metabolome_data) <- paste0("metabolite_", 1:p)
  return(list(metabolome_data = metabolome_data,
              metabolome_base = metabolome_base,
              micro_effect = micro_effect,
              metabolome_error = metabolome_error))
}

## example ##
effect_index(n = 100, q = 200, p = 200, seed = 1)
microbiome_scale <- apply(microbiome_relative_abundance, 2, X_scale_0)
metabolome_data <- metabolome_data_generate(n = 100,  n_effect = 100, q = 200, p = 200, effect_scale = 0.5, 
                                            microbiome_matrix = microbiome_scale, seed = 1)$metabolome_data
