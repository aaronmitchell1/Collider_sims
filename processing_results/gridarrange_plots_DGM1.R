library(tidyverse)
library(gridExtra)

#Plots for DGM 1 with adjusted y-axis limits

plot_mean_bias <- ggplot(summary_df, aes(x = Scenario, y = Mean_Bias, group = Method, color = Method, linetype = Method)) +
    geom_line(size = 1) + 
    geom_point(size = 3) + 
    labs(
        x = "Scenario",
        y = "Mean Bias") +
    scale_y_continuous(limits = c(-0.04, 0.01)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(drop = FALSE)

plot_mean_se <- ggplot(summary_df, aes(x = Scenario, y = Mean_SE, group = Method, color = Method, linetype = Method)) +
    geom_line(size = 1) + 
    geom_point(size = 3) + 
    labs(
        x = "Scenario",
        y = "Mean SE") +
    scale_y_continuous(limits = c(0, 0.1)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(drop = FALSE)

plot_empirical_coverage <- ggplot(summary_df, aes(x = Scenario, y = Empirical_Coverage, group = Method, color = Method, linetype = Method)) +
    geom_line(size = 1) + 
    geom_point(size = 3) + 
    labs(
        x = "Scenario",
        y = "Empirical Coverage") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(drop = FALSE)

plot_type_i_error_null <- ggplot(summary_df, aes(x = Scenario, y = Type_I_Error_Null, group = Method, color = Method, linetype = Method)) +
    geom_line(size = 1) + 
    geom_point(size = 3) + 
    labs(
        x = "Scenario",
        y = "T1E Rate") +
    theme_minimal() +
    scale_y_continuous(limits = c(0, 0.6)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(drop = FALSE)

plot_power_to_detect_nonnull <- ggplot(summary_df, aes(x = Scenario, y = Power_to_Detect_NonNull, group = Method, color = Method, linetype = Method)) +
    geom_line(size = 1) + 
    geom_point(size = 3) + 
    labs(
        x = "Scenario",
        y = "Power") +
    theme_minimal() +
    scale_y_continuous(limits = c(0.5, 1.0)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(drop = FALSE)

#Combine plots into one
combined_plot <- grid.arrange(
    plot_mean_bias + theme(legend.position = "none"), 
    plot_mean_se + theme(legend.position = "none"),
    plot_empirical_coverage + theme(legend.position = "none"),
    plot_type_i_error_null + theme(legend.position = "none"),
    plot_power_to_detect_nonnull + theme(legend.position = "none"),
    ncol = 2,
    top = ""  #This needs to be included and blank or gives an error.
)

#Add shared legend
legend <- get_legend(plot_mean_bias)  #Extract legend from any plot

#Combine plots and legend
plot <- grid.arrange(
    combined_plot,
    legend,
    ncol = 2,
    widths = c(4, 1)  #Adjust widths to add legend at the side
)

ggsave(plot, filename = "DGM1.png")
