summary_df <- data.frame(
  Method = c("Unadjusted", "Dudbridge", "Median", "MR-Raps", "MR-Horse", "SlopeHunter"),
  Mean_Bias = c(mean(true_mat - progression_mat), mean(true_mat - dudbridge_mat), mean(true_mat - median_mat), mean(true_mat - mr_raps_mat), mean(true_mat - mr_horse_mat), mean(true_mat - slopehunter_mat)),
  Mean_SE = c(mean(progression_se_mat), mean(dudbridge_se_mat), mean(median_se_mat), mean(mr_raps_se_mat), mean(mr_horse_se_mat), mean(slopehunter_se_mat)),
  Empirical_Coverage = c(emp.coverage(progression_mat, progression_se_mat, true_mat), emp.coverage(dudbridge_mat, dudbridge_se_mat, true_mat), emp.coverage(median_mat, median_se_mat, true_mat), emp.coverage(mr_raps_mat, mr_raps_se_mat, true_mat), emp.coverage(mr_horse_mat, mr_horse_se_mat, true_mat), emp.coverage(slopehunter_mat, slopehunter_se_mat, true_mat)),
  Type_I_Error_Null = c(1 - emp.coverage(progression_mat[, 1:80], progression_se_mat[, 1:80], 0), 1 - emp.coverage(dudbridge_mat[, 1:80], dudbridge_se_mat[, 1:80], 0), 1 - emp.coverage(median_mat[, 1:80], median_se_mat[, 1:80], 0), 1 - emp.coverage(mr_raps_mat[, 1:80], mr_raps_se_mat[, 1:80], 0), 1 - emp.coverage(mr_horse_mat[, 1:80], mr_horse_se_mat[, 1:80], 0), 1 - emp.coverage(slopehunter_mat[, 1:80], slopehunter_se_mat[, 1:80], 0)),
  Power_to_Detect_NonNull = c(1 - emp.coverage(progression_mat[, 81:100], progression_se_mat[, 81:100], 0), 1 - emp.coverage(dudbridge_mat[, 81:100], dudbridge_se_mat[, 81:100], 0), 1 - emp.coverage(median_mat[, 81:100], median_se_mat[, 81:100], 0), 1 - emp.coverage(mr_raps_mat[, 81:100], mr_raps_se_mat[, 81:100], 0), 1 - emp.coverage(mr_horse_mat[, 81:100], mr_horse_se_mat[, 81:100], 0), 1 - emp.coverage(slopehunter_mat[, 81:100], slopehunter_se_mat[, 81:100], 0))
)

p1 <- ggplot(summary_df, aes(x = Method, y = Mean_Bias)) +
geom_bar(stat = "identity", fill = "skyblue") +
theme_minimal() +
labs(y = "Mean Bias", x = "") +
theme(axis.text.x = element_text(angle = 45, hjust = 1))

p2 <- ggplot(summary_df, aes(x = Method, y = Mean_SE)) +
geom_bar(stat = "identity", fill = "salmon") +
theme_minimal() +
labs(y = "Mean SE", x = "") +
scale_y_continuous(limits = c(0, 0.03)) +
theme(axis.text.x = element_text(angle = 45, hjust = 1))

p3 <- ggplot(summary_df, aes(x = Method, y = Empirical_Coverage)) +
geom_bar(stat = "identity", fill = "lightgreen") +
theme_minimal() +
labs(y = "Coverage", x = "") +
theme(axis.text.x = element_text(angle = 45, hjust = 1))

p4 <- ggplot(summary_df, aes(x = Method, y = Type_I_Error_Null)) +
geom_bar(stat = "identity", fill = "navyblue") +
theme_minimal() +
labs(y = "T1E Rate", x = "") +
scale_y_continuous(limits = c(0, 0.4)) +
theme(axis.text.x = element_text(angle = 45, hjust = 1))

p5 <- ggplot(summary_df, aes(x = Method, y = Power_to_Detect_NonNull)) +
geom_bar(stat = "identity", fill = "lightblue") +
theme_minimal() +
labs(y = "Power", x = "") +
scale_y_continuous(limits = c(0, 1.0)) +
theme(axis.text.x = element_text(angle = 45, hjust = 1))

g1 <- grid.arrange(p1, p2, p3, p4, p5, ncol = 2)
ggsave(file="dgm1.png", g1)
