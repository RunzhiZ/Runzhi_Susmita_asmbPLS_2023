######## PARAMETERS ########**********************************************
#### n: number of samples
#### q: number of microbial features
#### m: number of baseline sequence reads for microbial data
#### seed: seed for data generation
####**********************************************************************

microbiome_data_generate <- function(n, q, m, seed) {
  set.seed(seed)
  ## value of alpha for dirichlet distribution,  here chosen at random
  alpha_dirichlet <- runif(q)
  microbiome_relative_abundance <- microbiome_raw_data <- matrix(NA, nrow = n,
                                                                 ncol = q)
  for (i in 1:n) {
    set.seed(i)
    ## library size
    l <- round(runif(1, m, 2 * m))
    ## pre-simulation of the Dirichlet
    proportion <- rgamma(q, alpha_dirichlet, 1)
    microbiome_raw_data[i, ] <- rmultinom(1, l, proportion)
  }
  microbiome_relative_abundance <- t(apply(microbiome_raw_data, 1, function(data) {
    data / sum(data)
  }))
  colnames(microbiome_raw_data) <- colnames(microbiome_relative_abundance) <- paste0("Taxa_", 1:q)
  return(list(microbiome_raw_data = microbiome_raw_data,
              microbiome_relative_abundance = microbiome_relative_abundance))
}

## example ##
system.time(microbiome_raw_data <- microbiome_data_generate(n = 100, q = 200, m = 2000, seed = 1)$microbiome_raw_data)
system.time(microbiome_relative_abundance <- microbiome_data_generate(n = 100, q = 200, m = 2000, seed = 1)$microbiome_relative_abundance)

