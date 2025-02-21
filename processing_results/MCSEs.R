#Vectors to store MCSE for each method
mcse_dudbridge <- numeric(reps)
mcse_weighted_median <- numeric(reps)
mcse_mr_raps <- numeric(reps)
mcse_mr_horse <- numeric(reps)
mcse_slopehunter <- numeric(reps)

#Loop through iterations, extract and store correction betas
for (n in 1:reps) {
  collider_bias_results <- loop_methods[[n]]$collider_bias_results
  
  dudbridge_betas <- collider_bias_results[collider_bias_results$Method == 'Dudbridge', 'Correction_Beta']
  weighted_median_betas <- collider_bias_results[collider_bias_results$Method == 'Weighted_median', 'Correction_Beta']
  mr_raps_betas <- collider_bias_results[collider_bias_results$Method == 'MR_RAPS', 'Correction_Beta']
  mr_horse_betas <- collider_bias_results[collider_bias_results$Method == 'MR_Horse', 'Correction_Beta']
  slopehunter_betas <- collider_bias_results[collider_bias_results$Method == 'Slopehunter', 'Correction_Beta']
  
  mcse_dudbridge[n] <- dudbridge_betas
  mcse_weighted_median[n] <- weighted_median_betas
  mcse_mr_raps[n] <- mr_raps_betas
  mcse_mr_horse[n] <- mr_horse_betas
  mcse_slopehunter[n] <- slopehunter_betas
}

#Calculate Monte Carlo standard errors for each method
mcse_dudbridge <- sd(mcse_dudbridge) / sqrt(reps)
mcse_weighted_median <- sd(mcse_weighted_median) / sqrt(reps)
mcse_mr_raps <- sd(mcse_mr_raps) / sqrt(reps)
mcse_mr_horse <- sd(mcse_mr_horse) / sqrt(reps)
mcse_slopehunter <- sd(mcse_slopehunter) / sqrt(reps)
