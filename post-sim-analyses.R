#Create a df from the results loop for easier analysis of performance metrics.

library(dplyr)
library(rsimsum)

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


#Set up data frame to store the true values of θ.

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

#Store incidence betas in a data frame.

incidence_df <- data.frame()

for (i in 1:reps) {
  incidence_results <- loop_incidence_GWAS[[i]]$incidence_GWAS
  temp_df <- tibble(
    Beta = as.numeric(incidence_results[, 1]),
    Iteration = i
  )
  incidence_df <- bind_rows(incidence_df, temp_df)
  
}

#Take the rolling mean of θ for nSNPs for each iteration.
mean_values <- as.data.frame(tapply(true_df$true_value, true_df$Iteration, mean))

#Mean beta for incidence SNPs for each iteration.
mean_betas <- as.data.frame(tapply(incidence_df$Beta, incidence_df$Iteration, mean))

#This is calculating bias in terms of θ for variables that cause selection.
#Generate the true values of θ for each method (unusually for a simulation this relies on knowing the value of the estimate, the correction factor slope) using Apostolis' formula:
#true β P (θ) = mean(var_interaction) + correction factor slope for each method * β I. 
#This method makes the assumptions that the model for collider is perfectly log-additive & exposure is binary.
#Will turn this into a function if more methods added.

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

theta_results <- cbind(results_df, true_theta_Dudbridge, true_theta_MR_RAPS, true_theta_Weighted_median, true_theta_MR_Horse, true_theta_SlopeHunter)

#The MR-Horse method is Bayesian so outputs SDs, to turn them into SEs to make them the same as other estimators I assume you use SE = SD/√(n individuals)
results_df$Correction_SE[results_df$Method=="MR_Horse"] <- (results_df$Correction_SE[results_df$Method=="MR_Horse"]/sqrt(n.ind))

#Put the methods into their own data frame to conduct the rsimsum analysis, I don't think you can yet have multiple methods with their own true value of θ in a single call to rsimusum.
results_data_Dudbridge <- data.frame(results_df$Correction_Beta[results_data$Method=="Dudbridge"], results_df$Correction_SE[results_data$Method=="Dudbridge"], results_df$Iteration[results_data$Method=="Dudbridge"], true_theta_Dudbridge); colnames(results_data_Dudbridge) = c("Estimate", "StdErr", "Iter", "True"))
results_data_Weighted_median <- data.frame(results_df$Correction_Beta[results_data$Method=="Weighted_median"], results_df$Correction_SE[results_data$Method=="Weighted_median"], results_df$Iteration[results_data$Method=="Weighted_median"], true_theta_Weighted_median); colnames(results_data_Weighted_median) = c("Estimate", "StdErr", "Iter", "True")
results_data_MR_RAPS <- data.frame(results_df$Correction_Beta[results_data$Method=="MR_RAPS"], results_df$Correction_SE[results_data$Method=="MR_RAPS"], results_df$Iteration[results_data$Method=="MR_RAPS"], true_theta_MR_RAPS); colnames(results_data_MR_RAPS) = c("Estimate", "StdErr", "Iter", "True")
results_data_MR_Horse <- data.frame(results_df$Correction_Beta[results_data$Method=="MR_Horse"], results_df$Correction_SE[results_data$Method=="MR_Horse"], results_df$Iteration[results_data$Method=="MR_Horse"], true_theta_MR_Horse); colnames(results_data_MR_Horse) = c("Estimate", "StdErr", "Iter", "True")
results_data_Slopehunter <- data.frame(results_df$Correction_Beta[results_data$Method=="Slopehunter"], results_df$Correction_SE[results_data$Method=="Slopehunter"], results_df$Iteration[results_data$Method=="Slopehunter"], true_theta_SlopeHunter); colnames(results_data_Slopehunter) = c("Estimate", "StdErr", "Iter", "True")

#Run rsimsum to get the performance metrics for each method.
res_Dudbridge <- simsum(data = results_data_Dudbridge, estvarname = "Estimate", se = "StdErr", true = "True")
res_Weighted_median <- simsum(data = results_data_Weighted_median, estvarname = "Estimate", se = "StdErr", true = "True")
res_MR_RAPS <- simsum(data = results_data_MR_RAPS, estvarname = "Estimate", se = "StdErr", true = "True")
res_MR_Horse <- simsum(data = results_data_MR_Horse, estvarname = "Estimate", se = "StdErr", true = "True")
res_Slopehunter <- simsum(data = results_data_Slopehunter, estvarname = "Estimate", se = "StdErr", true = "True")
