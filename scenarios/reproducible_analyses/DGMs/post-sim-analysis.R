library(dplyr)
library(rsimsum)
reps <- 1000

#Get methods results in a dataframe.

results_df <- data.frame()

for (n in 1:reps) {
  # Extract results for each method for the current iteration
  collider_results <- loop_methods[[n]]$collider_bias_results
  
  # Create a temporary dataframe for the current iteration results
  temp_df <- tibble(
    Method = as.character(collider_results$Method),
    Correction_Beta = as.numeric(collider_results$Correction_Beta),
    Correction_SE = as.numeric(collider_results$Correction_SE),
    Iteration = n
  )
  
  # Bind the temporary dataframe to the overall results dataframe
  results_df <- bind_rows(results_df, temp_df)
  
}

#Get true values in a dataframe.

true_df <- data.frame()

for (n in 1:reps) {
  # Extract results for each method for the current iteration
  true_results <- loop_var_interaction[[n]]$interaction_df...var.U.
  
  # Create a temporary dataframe for the current iteration results
  temp_df <- tibble(
    true_value = as.numeric(true_results),
    Iteration = n
  )
  
  # Bind the temporary dataframe to the overall results dataframe
  true_df <- bind_rows(true_df, temp_df)
  
}

#Get incidence betas in a dataframe.

incidence_df <- data.frame()

for (i in 1:reps) {
  incidence_results <- loop_incidence_GWAS[[i]]$incidence_GWAS
  temp_df <- tibble(
    Beta = as.numeric(incidence_results[, 1]),
    Iteration = i
  )
  incidence_df <- bind_rows(incidence_df, temp_df)
  
}

#Get progression betas in a dataframe.

progression_df <- data.frame()

for (i in 1:reps) {
  progression_results <- loop_progression_GWAS[[i]]$progression_GWAS
  temp_df <- tibble(
    Beta = as.numeric(progression_results[, 1]),
    Iteration = i
  )
  progression_df <- bind_rows(progression_df, temp_df)
  
}

pobs_df <- data.frame()

for (i in 1:reps) {
  pobs_results <- loop_progression_oracle_GWAS[[i]]$progression_oracle_GWAS
  temp_df <- tibble(
    Beta = as.numeric(pobs_results[, 1]),
    Iteration = i
  )
  pobs_df <- bind_rows(pobs_df, temp_df)
  
}

#Take the mean betas for P and Pobs across all reps and get the difference.

mean_p_betas <- as.data.frame(tapply(progression_df$Beta, progression_df$Iteration, mean)); colnames(mean_p_betas) <- "Mean"
mean_pobs_betas <- as.data.frame(tapply(pobs_df$Beta, pobs_df$Iteration, mean)); colnames(mean_pobs_betas) <- "Mean"
mean(mean_pobs_betas$Mean) - mean(mean_p_betas$Mean)

#Take the means of theta for nSNPs for each iteration
mean_values <- as.data.frame(tapply(true_df$true_value, true_df$Iteration, mean))

#Mean beta for incidence SNPs for each iteration
mean_betas <- as.data.frame(tapply(incidence_df$Beta, incidence_df$Iteration, mean))

#Apostolis's formula for true value for each method: this is calculating bias in terms of theta for variables that cause selection. Assumptions:
#model for collider is perfectly log-additive, exposure binary, 

true_theta_Dudbridge <- data.frame()

true_theta_Dudbridge <- as.data.frame(mean_values$`tapply(true_df$true_value, true_df$Iteration, mean)` + results_df$Correction_Beta[results_df$Method=="Dudbridge"] * mean_betas)

true_theta_Weighted_median <- data.frame()
  
true_theta_Weighted_median <- as.data.frame(mean_values$`tapply(true_df$true_value, true_df$Iteration, mean)` + results_df$Correction_Beta[results_df$Method=="Weighted_median"] * mean_betas)

true_theta_MR_RAPS <- data.frame()
    
true_theta_MR_RAPS <- as.data.frame(mean_values$`tapply(true_df$true_value, true_df$Iteration, mean)` + results_df$Correction_Beta[results_df$Method=="MR_RAPS"] * mean_betas)

true_theta_MR_Horse <- data.frame()

true_theta_MR_Horse <- as.data.frame(mean_values$`tapply(true_df$true_value, true_df$Iteration, mean)` + results_df$Correction_Beta[results_df$Method=="MR_Horse"] * mean_betas)

true_theta_SlopeHunter <- data.frame()

true_theta_SlopeHunter <- as.data.frame(mean_values$`tapply(true_df$true_value, true_df$Iteration, mean)` + results_df$Correction_Beta[results_df$Method=="Slopehunter"] * mean_betas)

results_data <- cbind(results_df, true_theta_Dudbridge, true_theta_MR_RAPS, true_theta_Weighted_median, true_theta_MR_Horse, true_theta_SlopeHunter)

results_data_Dudbridge <- data.frame(results_df$Correction_Beta[results_data$Method=="Dudbridge"], results_df$Correction_SE[results_data$Method=="Dudbridge"], results_df$Iteration[results_data$Method=="Dudbridge"], true_theta_Dudbridge); colnames(results_data_Dudbridge) = c("Estimate", "StdErr", "Iter", "True")
results_data_Weighted_median <- data.frame(results_df$Correction_Beta[results_data$Method=="Weighted_median"], results_df$Correction_SE[results_data$Method=="Weighted_median"], results_df$Iteration[results_data$Method=="Weighted_median"], true_theta_Weighted_median); colnames(results_data_Weighted_median) = c("Estimate", "StdErr", "Iter", "True")
results_data_MR_RAPS <- data.frame(results_df$Correction_Beta[results_data$Method=="MR_RAPS"], results_df$Correction_SE[results_data$Method=="MR_RAPS"], results_df$Iteration[results_data$Method=="MR_RAPS"], true_theta_MR_RAPS); colnames(results_data_MR_RAPS) = c("Estimate", "StdErr", "Iter", "True")
results_data_MR_Horse <- data.frame(results_df$Correction_Beta[results_data$Method=="MR_Horse"], results_df$Correction_SE[results_data$Method=="MR_Horse"], results_df$Iteration[results_data$Method=="MR_Horse"], true_theta_MR_Horse); colnames(results_data_MR_Horse) = c("Estimate", "StdErr", "Iter", "True")
results_data_Slopehunter <- data.frame(results_df$Correction_Beta[results_data$Method=="Slopehunter"], results_df$Correction_SE[results_data$Method=="Slopehunter"], results_df$Iteration[results_data$Method=="Slopehunter"], true_theta_SlopeHunter); colnames(results_data_Slopehunter) = c("Estimate", "StdErr", "Iter", "True")

res_Dudbridge <- simsum(data = results_data_Dudbridge, estvarname = "Estimate", se = "StdErr", true = "True")
res_Weighted_median <- simsum(data = results_data_Weighted_median, estvarname = "Estimate", se = "StdErr", true = "True")
res_MR_RAPS <- simsum(data = results_data_MR_RAPS, estvarname = "Estimate", se = "StdErr", true = "True")
res_MR_Horse <- simsum(data = results_data_MR_Horse, estvarname = "Estimate", se = "StdErr", true = "True")
res_Slopehunter <- simsum(data = results_data_Slopehunter, estvarname = "Estimate", se = "StdErr", true = "True")

#Get the difference between methods betas and P complete cases betas.

db <- subset(results_df, Method == "Dudbridge")
med <- subset(results_df, Method == "Weighted_median")
RAPS <- subset(results_df, Method == "MR-RAPS")
horse <- subset(results_df, Method == "MR_Horse")
slopehunter <- subset(results_df, Method == "Slopehunter")

mean(db$Correction_Beta) - mean(mean_p_betas$Mean)
mean(med$Correction_Beta) - mean(mean_p_betas$Mean)
mean(RAPS$Correction_Beta) - mean(mean_p_betas$Mean)
mean(horse$Correction_Beta) - mean(mean_p_betas$Mean)
mean(slopehunter$Correction_Beta) - mean(mean_p_betas$Mean)

#Select only progression SNPs for b methods - P analysis (subtracting mean P complete cases betas from methods results on Pobs),
#N.B. the df xx:xx value needs to change depending on allocation of SNPs to clusters specified in DGM.

selected_rows <- list()
for(i in 1:reps) {df <- loop_progression_oracle_GWAS[[i]]$progression_oracle_GWAS
   selected <- df[81:90, 1]
   selected_rows[[i]] <- selected
 }

#Put all progression only SNPs from all reps into a new dataframe.

prog_SNPs_only <- do.call(rbind, selected_rows)

#Calculate b methods - P complete cases for only progression SNPs.

mean(db$Correction_Beta) - mean(prog_SNPs_only)
mean(med$Correction_Beta) - mean(prog_SNPs_only)
mean(RAPS$Correction_Beta) - mean(prog_SNPs_only)
mean(horse$Correction_Beta) - mean(prog_SNPs_only)
mean(slopehunter$Correction_Beta) - mean(prog_SNPs_only)
