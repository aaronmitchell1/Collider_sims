med_plot <- results_df %>%
    filter(Method == "Weighted_median") %>%
    ggplot(aes(x = Correction_Beta)) +
    geom_histogram(bins = 30, fill = "blue", color = "black", alpha = 0.7) +
    theme_minimal() +
    xlim(-10, 10) +
    labs(
        x = "Estimate of correction factor beta",
        y = "Reps"
    )

raps_plot <- results_df %>%
    filter(Method == "MR_RAPS") %>% 
    ggplot(aes(x = Correction_Beta)) + 
    geom_histogram(bins = 30, fill = "red", color = "black", alpha = 0.7) + 
    theme_minimal() + 
    ylim(0, 1000) +
    xlim(-10, 10) +
    labs(
        x = "Estimate of correction factor beta",
        y = "Reps"
    )

p1 <- grid.arrange(raps_plot, med_plot)
ggsave(p1, file = "rapsmed.png")
