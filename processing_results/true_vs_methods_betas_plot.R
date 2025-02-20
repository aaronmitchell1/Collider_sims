#Convert matrices back to long dfs for ggplot
true_df_long <- as.data.frame(true_mat) %>%
    mutate(Iteration = 1:reps) %>%
    gather(key = "SNP", value = "True_Beta", -Iteration)

dudbridge_df_long <- as.data.frame(dudbridge_mat) %>%
    mutate(Iteration = 1:reps) %>%
    gather(key = "SNP", value = "Dudbridge_Beta", -Iteration)

horse_df_long <- as.data.frame(mr_horse_mat) %>%
    mutate(Iteration = 1:reps) %>%
    gather(key = "SNP", value = "Horse_Beta", -Iteration)

median_df_long <- as.data.frame(median_mat) %>%
    mutate(Iteration = 1:reps) %>%
    gather(key = "SNP", value = "Median_Beta", -Iteration)

raps_df_long <- as.data.frame(mr_raps_mat) %>%
    mutate(Iteration = 1:reps) %>%
    gather(key = "SNP", value = "Raps_Beta", -Iteration)

slopehunter_df_long <- as.data.frame(slopehunter_mat) %>%
    mutate(Iteration = 1:reps) %>%
    gather(key = "SNP", value = "SlopeHunter_Beta", -Iteration)

#Merge the true data with each method's dfs by rep and SNP
merged_dudbridge <- merge(true_df_long, dudbridge_df_long, by = c("Iteration", "SNP"))
merged_horse <- merge(true_df_long, horse_df_long, by = c("Iteration", "SNP"))
merged_median <- merge(true_df_long, median_df_long, by = c("Iteration", "SNP"))
merged_raps <- merge(true_df_long, raps_df_long, by = c("Iteration", "SNP"))
merged_slopehunter <- merge(true_df_long, slopehunter_df_long, by = c("Iteration", "SNP"))

plot_dudbridge <- ggplot(merged_dudbridge, aes(x = True_Beta, y = Dudbridge_Beta)) +
    geom_point(alpha = 0.5) + 
    geom_abline(slope = 1, intercept = 0, color = "red") +
    theme_minimal() +
    labs(x = "True beta", y = "Dudbridge beta")

plot_horse <- ggplot(merged_horse, aes(x = True_Beta, y = Horse_Beta)) +
    geom_point(alpha = 0.5) + 
    geom_abline(slope = 1, intercept = 0, color = "red") +
    theme_minimal() +
    labs(x = "True beta", y = "MR-Horse beta")

plot_median <- ggplot(merged_median, aes(x = True_Beta, y = Median_Beta)) +
    geom_point(alpha = 0.5) + 
    geom_abline(slope = 1, intercept = 0, color = "red") +
    theme_minimal() +
    labs(x = "True beta", y = "Weighted median beta")

plot_raps <- ggplot(merged_raps, aes(x = True_Beta, y = Raps_Beta)) +
    geom_point(alpha = 0.5) + 
    geom_abline(slope = 1, intercept = 0, color = "red") +
    theme_minimal() +
    labs(x = "True beta", y = "MR-RAPS beta") +
    scale_y_continuous(limits = c(-1.5, 1.5)) +
    scale_x_continuous(limits = c(-1.5, 1.5))

plot_slopehunter <- ggplot(merged_slopehunter, aes(x = True_Beta, y = SlopeHunter_Beta)) +
    geom_point(alpha = 0.5) + 
    geom_abline(slope = 1, intercept = 0, color = "red") +
    theme_minimal() +
    labs(x = "True beta", y = "SlopeHunter beta")

p1 <- grid.arrange(plot_dudbridge, plot_horse, plot_median, plot_raps, plot_slopehunter, ncol = 2)
ggsave(p1, file = "truevsmethodsbetas.png")
