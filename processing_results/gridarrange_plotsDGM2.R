library(patchwork)

p1 <- ggplot(summary_df, aes(x = Scenario, y = Mean_Bias, color = Method, group = Method)) +
geom_line(size = 1) +  
geom_point(size = 3) +
theme_minimal() +    
labs(x = "Scenario",
y = "Mean Bias") +
coord_cartesian(ylim=c(-0.00005, 0.0002)) +
theme(axis.text.x = element_text(angle = 45, hjust = 1),
legend.position = "bottom")

p2 <- ggplot(summary_df, aes(x = Scenario, y = Mean_SE, color = Method, group = Method)) +
geom_line(size = 1) +
geom_point(size = 3) +
theme_minimal() +
coord_cartesian(ylim=c(0.0045, 0.006)) +
labs(x = "Scenario",
y = "Mean SE") +
theme(axis.text.x = element_text(angle = 45, hjust = 1),
legend.position = "bottom")

p3 <- ggplot(summary_df, aes(x = Scenario, y = Type_I_Error_Null, color = Method, group = Method)) +
geom_line(size = 1) + geom_point(size = 3) + theme_minimal() +
labs(x = "Scenario", y = "Type I Error") +
scale_y_continuous(limits = c(0, 0.4)) +  
theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")

p4 <- ggplot(summary_df, aes(x = Scenario, y = Power_to_Detect_NonNull, color = Method, group = Method)) +
geom_line(size = 1) + geom_point(size = 3) + theme_minimal() +
labs(x = "Scenario", y = "Power") +
scale_y_continuous(limits = c(0.8, 1.0)) +  
theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")

p5 <- ggplot(summary_df, aes(x = Scenario, y = Empirical_Coverage, color = Method, group = Method)) +
geom_line(size = 1) + geom_point(size = 3) + theme_minimal() +
labs(x = "Scenario", y = "Empirical Coverage") +
scale_y_continuous(limits = c(0.7, 1.0)) +  
theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")

combined_plot <- (p1 | p2) / (p3 | p4) / (p5)
combined_plot + plot_layout(guides = "collect") & theme(legend.position = "bottom")
p6 <- combined_plot + plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave(p6, filename = "DGM2.png")
