#Process the results for DGM 1 where nSNPs was varied and make graph for bias and average SE.

reps <- 1000
library(readr)
library(ggplot2)
library(tidyr)
library(gridExtra)

scenarios <- c("801010", "701020", "601030", "501040", "401050")
path <- "/Users/vc23656/Downloads/Collider_sims/"

scenario_results <- list()

#Read in data for each scenario
for (scenario in scenarios) {
    file_path <- paste0(path, scenario, "/sim_results.RData")
    load(file_path)
    
    results_df <- data.frame()
    incidence_df <- data.frame()
    progression_df <- data.frame()
    true_df <- data.frame()
    
    for (n in 1:reps) {
        collider_results <- loop_methods[[n]]$collider_bias_results
        temp_df <- tibble(
            Method = as.character(collider_results$Method),
            Correction_Beta = as.numeric(collider_results$Correction_Beta),
            Correction_SE = as.numeric(collider_results$Correction_SE),
            Iteration = n,
            Scenario = scenario
        )
        results_df <- bind_rows(results_df, temp_df)
    }
    
    for (i in 1:reps) {
        incidence_results <- loop_incidence_GWAS[[i]]$incidence_GWAS
        temp_df <- tibble(
            Beta = as.numeric(incidence_results[, 1]),
            SE = as.numeric(incidence_results[, 2]),
            Iteration = i,
            Scenario = scenario
        )
        incidence_df <- bind_rows(incidence_df, temp_df)
    }
    
    for (i in 1:reps) {
        progression_results <- loop_progression_GWAS[[i]]$progression_GWAS
        temp_df <- tibble(
            Beta = as.numeric(progression_results[, 1]),
            SE = as.numeric(progression_results[, 2]),
            Iteration = i,
            Scenario = scenario
        )
        progression_df <- bind_rows(progression_df, temp_df)
    }
    
    for (i in 1:reps) {
        true_results <- loop_progression_oracle_GWAS[[i]]$progression_oracle_GWAS
        temp_df <- tibble(
            Beta = as.numeric(true_results[, 1]),
            SE = as.numeric(true_results[, 2]),
            Iteration = i,
            Scenario = scenario
        )
        true_df <- bind_rows(true_df, temp_df)
    }
    
    #Save results
    scenario_results[[scenario]] <- list(
        results_df = results_df,
        incidence_df = incidence_df,
        progression_df = progression_df,
        true_df = true_df
    )
}

#Function to calculate adjusted summary statistics per method
compute_adjusted_stats <- function(progression_mat, incidence_mat, results_mat, progression_se_mat, incidence_se_mat) {
    dudbridge_mat <- progression_mat - matrix(results_mat[, 1], nrow = reps, ncol = nSNPs) * incidence_mat
    median_mat <- progression_mat - matrix(results_mat[, 2], nrow = reps, ncol = nSNPs) * incidence_mat
    mr_raps_mat <- progression_mat - matrix(results_mat[, 3], nrow = reps, ncol = nSNPs) * incidence_mat
    mr_horse_mat <- progression_mat - matrix(results_mat[, 4], nrow = reps, ncol = nSNPs) * incidence_mat
    slopehunter_mat <- progression_mat - matrix(results_mat[, 5], nrow = reps, ncol = nSNPs) * incidence_mat
    
    # Compute adjusted standard errors per SNP
    dudbridge_se_mat <- sqrt(progression_se_mat^2 + (results_mat[, 1]^2 * incidence_se_mat^2) + (incidence_mat^2 * progression_se_mat^2) + (incidence_se_mat^2 * progression_se_mat^2))
    median_se_mat <- sqrt(progression_se_mat^2 + (results_mat[, 2]^2 * incidence_se_mat^2) + (incidence_mat^2 * progression_se_mat^2) + (incidence_se_mat^2 * progression_se_mat^2))
    mr_raps_se_mat <- sqrt(progression_se_mat^2 + (results_mat[, 3]^2 * incidence_se_mat^2) + (incidence_mat^2 * progression_se_mat^2) + (incidence_se_mat^2 * progression_se_mat^2))
    mr_horse_se_mat <- sqrt(progression_se_mat^2 + (results_mat[, 4]^2 * incidence_se_mat^2) + (incidence_mat^2 * progression_se_mat^2) + (incidence_se_mat^2 * progression_se_mat^2))
    slopehunter_se_mat <- sqrt(progression_se_mat^2 + (results_mat[, 5]^2 * incidence_se_mat^2) + (incidence_mat^2 * progression_se_mat^2) + (incidence_se_mat^2 * progression_se_mat^2))
    
    return(list(
        dudbridge_mat = dudbridge_mat, 
        median_mat = median_mat, 
        mr_raps_mat = mr_raps_mat, 
        mr_horse_mat = mr_horse_mat, 
        slopehunter_mat = slopehunter_mat,
        dudbridge_se_mat = dudbridge_se_mat, 
        median_se_mat = median_se_mat,
        mr_raps_se_mat = mr_raps_se_mat,
        mr_horse_se_mat = mr_horse_se_mat, 
        slopehunter_se_mat = slopehunter_se_mat
    ))
}

#Calculate bias for each rep and method and summarise 

results_summary <- data.frame()

for (scenario in scenarios) {
    results_df <- scenario_results[[scenario]]$results_df
    incidence_df <- scenario_results[[scenario]]$incidence_df
    progression_df <- scenario_results[[scenario]]$progression_df
    true_df <- scenario_results[[scenario]]$true_df
    
    #Create matrices for reps, Betas and SEs
    progression_mat <- matrix(progression_df$Beta, nrow = reps, byrow = TRUE)
    incidence_mat <- matrix(incidence_df$Beta, nrow = reps, byrow = TRUE)
    true_mat <- matrix(true_df$Beta, nrow = reps, byrow = TRUE)
    
    progression_se_mat <- matrix(progression_df$SE, nrow = reps, byrow = TRUE)
    incidence_se_mat <- matrix(incidence_df$SE, nrow = reps, byrow = TRUE)
    true_se_mat <- matrix(true_df$SE, nrow = reps, byrow = TRUE)
    
    results_mat <- matrix(results_df$Correction_Beta, nrow = reps, byrow = TRUE)
    results_se_mat <- matrix(results_df$Correction_SE, nrow = reps, byrow = TRUE)
    
    #Calculate adjusted summary stats
    stats <- compute_adjusted_stats(
        progression_mat, 
        incidence_mat, 
        results_mat, 
        progression_se_mat, 
        incidence_se_mat
    )
    
    #Calculate bias (True - Estimated) for each rep and method
    for (method in c("Dudbridge", "Median", "MR-Raps", "MR-Horse", "SlopeHunter")) {
        method_col <- switch(method,
                             "Dudbridge" = stats$dudbridge_mat,
                             "Median" = stats$median_mat,
                             "MR-Raps" = stats$mr_raps_mat,
                             "MR-Horse" = stats$mr_horse_mat,
                             "SlopeHunter" = stats$slopehunter_mat
        )
        
        #Calculate bias (True - Estimated)
        bias <- true_mat - method_col
        
        #Average bias per SNP across all 1000 replications
        avg_bias <- rowMeans(bias)
        
        #Calculate SE for each method
        se <- switch(method,
                     "Dudbridge" = stats$dudbridge_se_mat,
                     "Median" = stats$median_se_mat,
                     "MR-Raps" = stats$mr_raps_se_mat,
                     "MR-Horse" = stats$mr_horse_se_mat,
                     "SlopeHunter" = stats$slopehunter_se_mat
        )
        
        #Average SE across all 1000 reps
        avg_se <- rowMeans(se)
        
        # Add to summary table for the current scenario
        results_summary <- bind_rows(results_summary, tibble(
            Scenario = scenario,
            Method = method,
            Avg_Bias = mean(avg_bias),  # Average bias across all SNPs
            Avg_SE = mean(avg_se)  # Average SE across all SNPs
        ))
    }
}

#Line graph for bias across reps.
p1 <- ggplot(results_summary, aes(x = Scenario, y = Avg_Bias, color = Method, group = Method)) +
    geom_line(linewidth = 1) +       #Use `linewidth` in newer version of ggplot
    geom_point(size = 3) +           
    theme_minimal() +
    labs(y = "Mean bias over 1000 replications", 
         x = "Scenario") +            
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),  #Rotate the X axis labels for readability
        legend.title = "Method",  
        plot.title = element_blank()     #No title
    )

#Line graph for average standard errors across reps.
p2 <- ggplot(results_summary, aes(x = Scenario, y = Avg_SE, color = Method, group = Method)) +
    geom_line(linewidth = 1) +       
    geom_point(size = 3) +           
    theme_minimal() +
    labs(y = "Mean Standard Error over 1000 replications", 
         x = "Scenario") +            
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.title = "Method",  
        plot.title = element_blank()     
    )

grid.arrange(grobs = list(p1, p2))
