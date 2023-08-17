######## PARAMETERS ########**********************************************
#### n: number of samples
####**********************************************************************

#### CENSORED IMPUTATION ####
#### Mean imputation
KM_estimator <- function(survival_data) {
  ## event indicator = 1 with no repeat
  survival_unique_order<-unique(survival_data[which(survival_data[,"Event_indicator"]==1), 1:2])[,"log_right_censored_time"]
  ## KM estimator of failure times
  KM_table <- matrix(NA, nrow = length(survival_unique_order), ncol = 5)
  colnames(KM_table) <- c("Time", "Number_of_failure", "Number_at_risk", "Survival_at_this_point", "Survival_total")
  for (i in 1:length(survival_unique_order)) {
    if(i == 1){
      previous_survival <- 1
    } else {
      previous_survival <- KM_table[i-1, "Survival_total"]
    }
    KM_table[i, "Time"] <- survival_unique_order[i]
    KM_table[i, "Number_of_failure"] <- length(which((survival_data[, "log_right_censored_time"] == KM_table[i, "Time"]) & 
                                                       (survival_data[, "Event_indicator"] == 1)))
    KM_table[i, "Number_at_risk"] <- length(which(survival_data[, "log_right_censored_time"] >= KM_table[i, "Time"]))
    KM_table[i, "Survival_at_this_point"] <- 1 - KM_table[i, "Number_of_failure"]/KM_table[i, "Number_at_risk"]
    KM_table[i, "Survival_total"] <- KM_table[i, "Survival_at_this_point"] * previous_survival
  }
  return(KM_table)
}

mean_imputation<-function(n, survival_data,round = F) {
  survival_data <- cbind(survival_data, 1:n, NA)
  colnames(survival_data)[3:4] <- c("ID", "Imputed_time")
  survival_data <- survival_data[order(survival_data[, "log_right_censored_time"]), ]
  survival_data[n, "Event_indicator"]<-1
  if (round == T) {
    survival_data[, "log_right_censored_time"] <- round(survival_data[, "log_right_censored_time"])
  }
  KM_table <- KM_estimator(survival_data)
  
  censored_index <- which(survival_data[, "Event_indicator"] == 0)
  
  for (i in 1:length(censored_index)) {
    censored_time <- survival_data[censored_index[i], "log_right_censored_time"]
    min_survival_time_index <- min(which(KM_table[, "Time"] > censored_time))
    numerator_table <- matrix(NA, nrow = nrow(KM_table) - min_survival_time_index + 1, ncol = 2)
    colnames(numerator_table) <- c("Tao", "Survival_diff")
    numerator_table[, "Tao"] <- KM_table[min_survival_time_index:nrow(KM_table), "Time"]
    if (min_survival_time_index == 1) {
      numerator_table[, "Survival_diff"] <- c(1, KM_table[min_survival_time_index:(nrow(KM_table)-1), "Survival_total"]) -
        KM_table[min_survival_time_index:nrow(KM_table), "Survival_total"]
    } else {
      numerator_table[, "Survival_diff"] <- KM_table[(min_survival_time_index-1):(nrow(KM_table)-1), "Survival_total"] -
        KM_table[min_survival_time_index:nrow(KM_table), "Survival_total"]
    }
    numerator <- sum(numerator_table[, "Tao"] * numerator_table[, "Survival_diff"])
    if (censored_time >= min(KM_table[, "Time"])) {
      denominator <- KM_table[max(which(KM_table[, "Time"] <= censored_time)), "Survival_total"]
    } else {
      denominator <- 1
    }
    survival_data[censored_index[i], "Imputed_time"] <- numerator / denominator
  }
  
  survival_data <- survival_data[order(survival_data[, "ID"]), ]
  survival_data[which(survival_data[,"Event_indicator"] == 1),"Imputed_time"] <- survival_data[which(survival_data[,"Event_indicator"] == 1),"log_right_censored_time"]
  return(list(imputed_table = survival_data,
              KM_table = KM_table))
}

## example ##
data_test <- matrix(c(1,1,1,2.5,5,7,1,1,0,1,0,1),ncol=2)
colnames(data_test) <- c("log_right_censored_time", "Event_indicator")
mean_imputation(6, data_test, F)
