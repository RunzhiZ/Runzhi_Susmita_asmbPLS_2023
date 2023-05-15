######## PARAMETERS ########**********************************************
#### n: number of samples
####**********************************************************************

#### CENSORED IMPUTATION ####
#### Reweighting
KM_estimator_censored <- function(survival_data) {
  ## event indicator = 1 with no repeat
  censored_unique_order<-unique(survival_data[which(survival_data[,"Event_indicator"]==0), 1:2])[,"log_right_censored_time"]
  ## KM estimator of censoring times
  KM_table_censored <- matrix(NA, nrow = length(censored_unique_order), ncol = 5)
  colnames(KM_table_censored) <- c("Time", "Number_censored", "Number_at_risk", "Censored_at_this_point", "Censored_total")
  for (i in 1:length(censored_unique_order)) {
    if(i == 1){
      previous_survival <- 1
    } else {
      previous_survival <- KM_table_censored[i-1, "Censored_total"]
    }
    KM_table_censored[i, "Time"] <- censored_unique_order[i]
    KM_table_censored[i, "Number_censored"] <- length(which((survival_data[, "log_right_censored_time"] == KM_table_censored[i, "Time"]) & 
                                                              (survival_data[, "Event_indicator"] == 0)))
    KM_table_censored[i, "Number_at_risk"] <- length(which(survival_data[, "log_right_censored_time"] >= KM_table_censored[i, "Time"])) - 
      length(which((survival_data[, "log_right_censored_time"] == KM_table_censored[i, "Time"]) & (survival_data[, "Event_indicator"] == 1)))
    KM_table_censored[i, "Censored_at_this_point"] <- 1 - KM_table_censored[i, "Number_censored"]/KM_table_censored[i, "Number_at_risk"]
    KM_table_censored[i, "Censored_total"] <- KM_table_censored[i, "Censored_at_this_point"] * previous_survival
  }
  return(KM_table_censored)
}

reweighting <- function(n, survival_data,round = F) {
  survival_data <- cbind(survival_data, 1:n, NA)
  colnames(survival_data)[3:4] <- c("ID", "Imputed_time")
  survival_data <- survival_data[order(survival_data[, "log_right_censored_time"]), ]
  survival_data[n, "Event_indicator"]<-1
  if (round == T) {
    survival_data[, "log_right_censored_time"] <- round(survival_data[, "log_right_censored_time"])
  }
  
  KM_table_censored <- KM_estimator_censored(survival_data)
  
  survival_index <- which(survival_data[, "Event_indicator"] == 1)
  
  for (i in which(survival_data[, "Event_indicator"] == 1)) {
    survival_time <- survival_data[i, "log_right_censored_time"]
    if (survival_time > min(KM_table_censored[, "Time"])) {
      survival_data[i, "Imputed_time"] <- survival_time / KM_table_censored[max(which(KM_table_censored[, "Time"] < survival_time)), "Censored_total"]
    } else {
      survival_data[i, "Imputed_time"] <- survival_time
    }
  }
  
  survival_data <- survival_data[order(survival_data[, "ID"]), ]
  survival_data[which(survival_data[, "Event_indicator"] == 0), "Imputed_time"] <- 0
  
  return(list(imputed_table = survival_data,
              KM_table_censored = KM_table_censored))
}

## example ##
data_test <- matrix(c(1,1,1,2.5,5,7,1,1,0,1,0,1),ncol=2)
colnames(data_test) <- c("log_right_censored_time", "Event_indicator")
reweighting(6, data_test, F)
sum(reweighting(n, survival_data[, 3:4], F)$imputed_table[, "Imputed_time"])/100
####******************************************************************###