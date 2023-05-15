setwd("C:/OneDrive/OneDrive - University of Florida/My research/Survival - microbiome/PLS/Data simulation")
load("Data Simulation function.RData")
load("Cencored data imputation function.RData")

#### parameters ####
n = 100
q = 150
p = 200
m = 1000
n_effect = 100
effect_scale = 0.5
r = 0.1
c= 0.3


#### Different strategies for beta ####
## First 5 - ,1:5,  remaining 0
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

#### data save ####
survival_data_1 <- one_step_simulation(n = 100, q = 150, p = 200, m = 2000,
                                       n_effect = 100, effect_scale = 0.5,
                                       r = 0.1, c = 0.3, beta = beta_1,
                                       seed = 2)
survival_data_2 <- one_step_simulation(n = 100, q = 150, p = 200, m = 2000,
                                       n_effect = 100, effect_scale = 0.5,
                                       r = 0.1, c = 0.3, beta = beta_2,
                                       seed = 2)
survival_data_3 <- one_step_simulation(n = 100, q = 150, p = 200, m = 2000,
                                       n_effect = 100, effect_scale = 0.5,
                                       r = 0.1, c = 0.3, beta = beta_3,
                                       seed = 2)
survival_data_4 <- one_step_simulation(n = 100, q = 150, p = 200, m = 2000,
                                       n_effect = 100, effect_scale = 0.5,
                                       r = 0.1, c = 0.3, beta = beta_4,
                                       seed = 2)

X <- survival_data_1$predictor_data

Y_1 <- matrix(mean_imputation(n, survival_data_1$survival_time_lognormal[, 3:4], F)$imputed_table[, 4], ncol = 1)
Y_2 <- matrix(mean_imputation(n, survival_data_2$survival_time_lognormal[, 3:4], F)$imputed_table[, 4], ncol = 1)
Y_3 <- matrix(mean_imputation(n, survival_data_3$survival_time_lognormal[, 3:4], F)$imputed_table[, 4], ncol = 1)
Y_4 <- matrix(mean_imputation(n, survival_data_4$survival_time_lognormal[, 3:4], F)$imputed_table[, 4], ncol = 1)

setwd("C:/OneDrive/OneDrive - University of Florida/My research/Survival - microbiome/PLS")
save(list=c('X',
            'Y_1',
            'Y_2',
            'Y_3',
            'Y_4'
),
file='Simulated_data.RData')
