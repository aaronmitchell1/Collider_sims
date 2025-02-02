#Varying the magnitude of index event bias (U*P effect) when generating SNP-progression associations.
#Format the average bias over 1000 reps for each method.

scenarios <- c("0.2", "0.4", "0.6", "0.8", "1.2", "1.4", "1.6", "1.8", "2.0")

#Empty dataframe to store results for all scenarios and methods
all_results <- data.frame()

#Loop through scenarios
for (scenario in scenarios) {
    
    file_path <- paste0("/Users/vc23656/Downloads/Collider_sims/", scenario, "/sim_results.RData")
    
    #Check if the file exists
    if (file.exists(file_path)) {
        #Load the simulation results for the current scenario
        load(file_path)
        
        #Get the number of replications from simulation results (1000)
        reps <- length(loop_methods)
        
        #Collect method results
        results_df <- data.frame()
        for (n in 1:reps) {
            collider_results <- loop_methods[[n]]$collider_bias_results
            temp_df <- tibble(
                Method = as.character(collider_results$Method),
                Correction_Beta = as.numeric(collider_results$Correction_Beta),
                Correction_SE = as.numeric(collider_results$Correction_SE),
                Iteration = n
            )
            results_df <- bind_rows(results_df, temp_df)
        }
        
        #Get incidence betas
        incidence_df <- data.frame()
        for (i in 1:reps) {
            incidence_results <- loop_incidence_GWAS[[i]]$incidence_GWAS
            temp_df <- tibble(
                Beta = as.numeric(incidence_results[, 1]),
                SE = as.numeric(incidence_results[, 2]),
                Iteration = i
            )
            incidence_df <- bind_rows(incidence_df, temp_df)
        }
        
        #Get progression betas
        progression_df <- data.frame()
        for (i in 1:reps) {
            progression_results <- loop_progression_GWAS[[i]]$progression_GWAS
            temp_df <- tibble(
                Beta = as.numeric(progression_results[, 1]),
                SE = as.numeric(progression_results[, 2]),
                Iteration = i
            )
            progression_df <- bind_rows(progression_df, temp_df)
        }
        
        #Get true progression betas in a dataframe
        true_df <- data.frame()
        for (i in 1:reps) {
            true_results <- loop_progression_oracle_GWAS[[i]]$progression_oracle_GWAS
            temp_df <- tibble(
                Beta = as.numeric(true_results[, 1]),
                SE = as.numeric(true_results[, 2]),
                Iteration = i
            )
            true_df <- bind_rows(true_df, temp_df)
        }
        
        #Create matrices for results
        progression_mat <- matrix(progression_df$Beta, nrow = reps, byrow = TRUE)
        incidence_mat <- matrix(incidence_df$Beta, nrow = reps, byrow = TRUE)
        true_mat <- matrix(true_df$Beta, nrow = reps, byrow = TRUE)
        
        #Correction matrices for methods
        results_mat <- matrix(results_df$Correction_Beta, nrow = reps, byrow = TRUE)
        results_se_mat <- matrix(results_df$Correction_SE, nrow = reps, byrow = TRUE)
        
        #Compute adjusted progression summary statistics for each method
        dudbridge_mat <- progression_mat - matrix(results_mat[, 1], reps, ncol(progression_mat)) * incidence_mat
        median_mat <- progression_mat - matrix(results_mat[, 2], reps, ncol(progression_mat)) * incidence_mat
        mr_raps_mat <- progression_mat - matrix(results_mat[, 3], reps, ncol(progression_mat)) * incidence_mat
        mr_horse_mat <- progression_mat - matrix(results_mat[, 4], reps, ncol(progression_mat)) * incidence_mat
        slopehunter_mat <- progression_mat - matrix(results_mat[, 5], reps, ncol(progression_mat)) * incidence_mat
        
        #Average bias calculation (difference between true and progression results)
        avg_bias <- data.frame(
            Scenario = scenario,
            Method = c("Dudbridge", "Median", "MR-Raps", "MR-Horse", "SlopeHunter"),
            Avg_Bias = c(mean(true_mat - dudbridge_mat),
                         mean(true_mat - median_mat),
                         mean(true_mat - mr_raps_mat),
                         mean(true_mat - mr_horse_mat),
                         mean(true_mat - slopehunter_mat)),
            Avg_SE = c(mean(results_se_mat[, 1]), mean(results_se_mat[, 2]), mean(results_se_mat[, 3]),
                       mean(results_se_mat[, 4]), mean(results_se_mat[, 5]))
        )
        
        #Combine results
        all_results <- bind_rows(all_results, avg_bias)
        
    } else {
        #If the file doesn't exist, print a warning and skip to the next scenario
        warning(paste("File not found for scenario", scenario, "at", file_path))
    }
}

#Make ggplot for bias.
ggplot(all_results, aes(x = Scenario, y = Avg_Bias, color = Method, group = Method)) +
    geom_line(linewidth = 1) +       # Use `linewidth` instead of `size`
    geom_point(size = 3) +           # Add points to the lines
    theme_minimal() +
    labs(y = "Mean bias over 1000 replications", 
         x = "Scenario") +            # Label for X axis
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate scenario names for readability
        legend.title = element_blank(),  # Optional: remove legend title
        plot.title = element_blank()     # Remove the title from the graph
    )

#Make ggplot for average standared errors.
ggplot(all_results, aes(x = Scenario, y = Avg_SE, color = Method, group = Method)) +
    geom_line(linewidth = 1) +       # Use `linewidth` instead of `size`
    geom_point(size = 3) +           # Add points to the lines
    theme_minimal() +
    labs(y = "Mean Standard Error over 1000 replications", 
         x = "Scenario") +            # Label for X axis
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate scenario names for readability
        legend.title = element_blank(),  # Optional: remove legend title
        plot.title = element_blank()     # Remove the title from the graph
    )
