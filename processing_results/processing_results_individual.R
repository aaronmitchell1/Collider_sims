#Set the scenario to "incb0.1" (90 I and 10 P SNPs, none affecting both)
scenario <- "incb0.1"

#Dataframe to store the results
summary_df <- tibble(
    Method = character(),
    Mean_Bias = numeric(),
    Mean_SE = numeric(),
    Empirical_Coverage = numeric(),
    Type_I_Error_Null = numeric(),
    Power_to_Detect_NonNull = numeric()
)

#Function to compute empirical coverage
emp.coverage <- function(est, se, true) {
    lb <- est - true + qnorm(0.025) * se
    ub <- est - true + qnorm(0.975) * se
    mean(lb * ub < 0)
}

#Load data for the specified scenario
load(paste0("/Users/vc23656/Downloads/Collider_sims/", scenario, "/sim_results.RData"))

#Convert lists to data frames for methods, incidence, progression, and true results
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

#Process incidence and progression results
incidence_df <- data.frame()
progression_df <- data.frame()
true_df <- data.frame()

for (i in 1:reps) {
    incidence_results <- loop_incidence_GWAS[[i]]$incidence_GWAS
    temp_df <- tibble(
        Beta = as.numeric(incidence_results[, 1]),
        SE = as.numeric(incidence_results[, 2]),
        Iteration = i
    )
    incidence_df <- bind_rows(incidence_df, temp_df)
    
    progression_results <- loop_progression_GWAS[[i]]$progression_GWAS
    temp_df <- tibble(
        Beta = as.numeric(progression_results[, 1]),
        SE = as.numeric(progression_results[, 2]),
        Iteration = i
    )
    progression_df <- bind_rows(progression_df, temp_df)
    
    true_results <- loop_progression_oracle_GWAS[[i]]$progression_oracle_GWAS
    temp_df <- tibble(
        Beta = as.numeric(true_results[, 1]),
        SE = as.numeric(true_results[, 2]),
        Iteration = i
    )
    true_df <- bind_rows(true_df, temp_df)
}

#Convert dataframes to matrices
progression_mat <- matrix(progression_df$Beta, nrow = reps, byrow = TRUE); colnames(progression_mat) <- paste("SNP", 1:nSNPs, sep = "")
incidence_mat <- matrix(incidence_df$Beta, nrow = reps, byrow = TRUE); colnames(incidence_mat) <- paste("SNP", 1:nSNPs, sep = "")
true_mat <- matrix(true_df$Beta, nrow = reps, byrow = TRUE); colnames(true_mat) <- paste("SNP", 1:nSNPs, sep = "")

#Standard error matrices
progression_se_mat <- matrix(progression_df$SE, nrow = reps, byrow = TRUE); colnames(progression_se_mat) <- paste("SNP", 1:nSNPs, sep = "")
incidence_se_mat <- matrix(incidence_df$SE, nrow = reps, byrow = TRUE); colnames(incidence_se_mat) <- paste("SNP", 1:nSNPs, sep = "")
true_se_mat <- matrix(true_df$SE, nrow = reps, byrow = TRUE); colnames(true_se_mat) <- paste("SNP", 1:nSNPs, sep = "")

#Results dataframe -> matrix
results_mat <- matrix(results_df$Correction_Beta, nrow = reps, byrow = TRUE)
colnames(results_mat) <- c("Dudbridge", "Median", "MR-Raps", "MR-Horse", "SlopeHunter")
results_se_mat <- matrix(results_df$Correction_SE, nrow = reps, byrow = TRUE)
colnames(results_se_mat) <- c("Dudbridge", "Median", "MR-Raps", "MR-Horse", "SlopeHunter")

#Adjusted matrices for each method
dudbridge_mat <- progression_mat - matrix(results_mat[, 1], reps, nSNPs) * incidence_mat
median_mat <- progression_mat - matrix(results_mat[, 2], reps, nSNPs) * incidence_mat
mr_raps_mat <- progression_mat - matrix(results_mat[, 3], reps, nSNPs) * incidence_mat
mr_horse_mat <- progression_mat - matrix(results_mat[, 4], reps, nSNPs) * incidence_mat
slopehunter_mat <- progression_mat - matrix(results_mat[, 5], reps, nSNPs) * incidence_mat

#Standard errors for adjusted matrices
dudbridge_se_mat <- sqrt(progression_se_mat^2 + (results_mat[, 1]^2 * incidence_se_mat^2) + (incidence_mat^2 * progression_se_mat^2) + (incidence_se_mat^2 * progression_se_mat^2))
median_se_mat <- sqrt(progression_se_mat^2 + (results_mat[, 2]^2 * incidence_se_mat^2) + (incidence_mat^2 * progression_se_mat^2) + (incidence_se_mat^2 * progression_se_mat^2))
mr_raps_se_mat <- sqrt(progression_se_mat^2 + (results_mat[, 3]^2 * incidence_se_mat^2) + (incidence_mat^2 * progression_se_mat^2) + (incidence_se_mat^2 * progression_se_mat^2))
mr_horse_se_mat <- sqrt(progression_se_mat^2 + (results_mat[, 4]^2 * incidence_se_mat^2) + (incidence_mat^2 * progression_se_mat^2) + (incidence_se_mat^2 * progression_se_mat^2))
slopehunter_se_mat <- sqrt(progression_se_mat^2 + (results_mat[, 5]^2 * incidence_se_mat^2) + (incidence_mat^2 * progression_se_mat^2) + (incidence_se_mat^2 * progression_se_mat^2))

#Calculate performance metrics for each method
method_results <- tibble(
    Method = c("Unadjusted", "Dudbridge", "Weighted Median", "MR-Raps", "MR-Horse", "SlopeHunter"),
    Mean_Bias = c(mean(true_mat - progression_mat), mean(true_mat - dudbridge_mat),
                  mean(true_mat - median_mat), mean(true_mat - mr_raps_mat),
                  mean(true_mat - mr_horse_mat), mean(true_mat - slopehunter_mat)),
    Mean_SE = c(mean(progression_se_mat), mean(dudbridge_se_mat), mean(median_se_mat),
                mean(mr_raps_se_mat), mean(mr_horse_se_mat), mean(slopehunter_se_mat)),
    Empirical_Coverage = c(
        emp.coverage(progression_mat, progression_se_mat, true_mat),
        emp.coverage(dudbridge_mat, dudbridge_se_mat, true_mat),
        emp.coverage(median_mat, median_se_mat, true_mat),
        emp.coverage(mr_raps_mat, mr_raps_se_mat, true_mat),
        emp.coverage(mr_horse_mat, mr_horse_se_mat, true_mat),
        emp.coverage(slopehunter_mat, slopehunter_se_mat, true_mat)
    ),
    Type_I_Error_Null = c(
        1 - emp.coverage(progression_mat[, 1:90], progression_se_mat[, 1:90], 0),
        1 - emp.coverage(dudbridge_mat[, 1:90], dudbridge_se_mat[, 1:90], 0),
        1 - emp.coverage(median_mat[, 1:90], median_se_mat[, 1:90], 0),
        1 - emp.coverage(mr_raps_mat[, 1:90], mr_raps_se_mat[, 1:90], 0),
        1 - emp.coverage(mr_horse_mat[, 1:90], mr_horse_se_mat[, 1:90], 0),
        1 - emp.coverage(slopehunter_mat[, 1:90], slopehunter_se_mat[, 1:90], 0)
    ),
    Power_to_Detect_NonNull = c(
        emp.coverage(progression_mat[, 91:nSNPs], progression_se_mat[, 91:nSNPs], true_mat[, 91:nSNPs]),
        emp.coverage(dudbridge_mat[, 91:nSNPs], dudbridge_se_mat[, 91:nSNPs], true_mat[, 91:nSNPs]),
        emp.coverage(median_mat[, 91:nSNPs], median_se_mat[, 91:nSNPs], true_mat[, 91:nSNPs]),
        emp.coverage(mr_raps_mat[, 91:nSNPs], mr_raps_se_mat[, 91:nSNPs], true_mat[, 91:nSNPs]),
        emp.coverage(mr_horse_mat[, 91:nSNPs], mr_horse_se_mat[, 91:nSNPs], true_mat[, 91:nSNPs]),
        emp.coverage(slopehunter_mat[, 91:nSNPs], slopehunter_se_mat[, 91:nSNPs], true_mat[, 91:nSNPs])
    )
)

#Store the results for the single scenario
summary_df <- bind_rows(summary_df, method_results)
